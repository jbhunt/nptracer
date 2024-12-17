import nrrd
import numpy as np
import pathlib as pl
from scipy import io
from scipy.spatial.transform import Rotation as R
from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.structure_tree import StructureTree

TREE = None
VOLUME = None
ANNOTATIONS = None
CCF_TO_STC = None
STC_TO_CCF = None

# TODO: Code this function
def _downloadCCF():
    """
    """

    return

def _loadAnnotations():
    """
    """

    packagePath = pl.Path(__file__).parent.parent
    ccfPath = packagePath.joinpath('resources', 'ccf')
    annotationsFile = ccfPath.joinpath('annotations.npy')
    if annotationsFile.exists() == False:
        raise Exception('Could not locate annotations file')
    global ANNOTATIONS
    ANNOTATIONS = np.load(annotationsFile)
    ANNOTATIONS = np.swapaxes(ANNOTATIONS, 1, 2)

    return

def _loadVolume():
    """
    """

    packagePath = pl.Path(__file__).parent.parent
    ccfPath = packagePath.joinpath('resources', 'ccf')
    volumeFile = ccfPath.joinpath('volume.npy')
    if volumeFile.exists() == False:
        raise Exception('Could not locate annotations file')
    global VOLUME
    VOLUME = np.load(volumeFile)
    VOLUME = np.swapaxes(VOLUME, 1, 2)

    return

def _loadAllenStructureTree():
    """
    Instatiate a brain structure tree object
    """

    oapi = OntologiesApi()
    structureGraph = oapi.get_structures_with_sets([1])  # 1 is the id of the adult mouse structure graph
    structureGraph = StructureTree.clean_structures(structureGraph)  
    global TREE
    TREE = StructureTree(structureGraph)

    return

def loadAllData():
    """
    """

    _loadVolume()
    _loadAnnotations()
    _loadAllenStructureTree()

    return

def translateBrainAreaIdentities(ids):
    """
    Identify the brain area acronym associated with each coded brain area
    identity
    """

    global TREE
    if TREE is None:
        _loadAllenStructureTree()
    structures = TREE.get_structures_by_id(ids)
    labels = [s['acronym'] if s is not None else '' for s in structures]

    return labels

def solveForElectrodeEndpoint(
    insertionPoint,
    insertionDepth=0,
    insertionAngle=15
    ):
    """
    Compute the endpoint of the electrode given the insertion point, the angle
    of insertion (0 to 180 degrees) and the insertion depth

    Notes
    -----
    The insertionAngle argument specifies the angle of insertion such that 0
    degrees means the probe is inserted purely left-to-right and 180 degrees
    means the probe is inserted purely right-to-left.
    """

    #
    B = np.copy(insertionPoint)

    try:
        assert insertionAngle >= 0
        assert insertionAngle <= 180
    except AssertionError:
        print('Insertion angle must be between 0 and 90 degrees')
        return

    #
    if insertionAngle == 0:
        B[1] -= insertionDepth

    elif insertionAngle == 90:
        B[2] += insertionDepth
 
    elif insertionAngle == 180:
        B[1] += insertionDepth

    elif insertionAngle > 90:
        insertionAngle = 180 - insertionAngle
        x = np.cos(np.deg2rad(insertionAngle)) * insertionDepth
        y = np.sin(np.deg2rad(insertionAngle)) * insertionDepth
        B[1] += x
        B[2] += y
        
    elif insertionAngle < 90:
        x = np.cos(np.deg2rad(insertionAngle)) * insertionDepth
        y = np.sin(np.deg2rad(insertionAngle)) * insertionDepth
        B[1] -= x
        B[2] += y

    return np.around(B, 2)

def getChannelPositionsForDenseConfig(
    dy=20,
    nChannels=374,
    sequence=[43, 11, 59, 27],
    ):
    """
    Create an array (N channels x 2) of the position of each electrode contact
    relative to the tip of the electrode (for Neuropixels 1.0 probes)
    """

    x = np.full(nChannels, np.nan)
    for i in range(nChannels):
        x[i] = np.take(sequence, i, mode='wrap')
    y = np.full(nChannels, np.nan)
    yi = dy
    for i in range(nChannels):
        if i % 2 != 0:
            yi += dy
        y[i] = yi

    xy = np.stack([x, y]).T

    return xy

def estimateSpikeDepths(workingDirectory):
    """
    Estimate the depth (distance from tip of electrode) of each spike
    """

    #
    if type(workingDirectory) == str:
        workingDirectory = pl.Path(workingDirectory)

    # Load data
    templateIndicesFile = workingDirectory.joinpath('spike_templates.npy')
    if templateIndicesFile.exists():
        templateIndices = np.load(templateIndicesFile)
    templatesFile = workingDirectory.joinpath('templates.npy')
    if templatesFile.exists():
        templates = np.load(templatesFile)
    channelPositionsFile = workingDirectory.joinpath('channel_positions.npy')
    if channelPositionsFile.exists():
        channelDepths = np.load(channelPositionsFile)[:, 1]

    # Compute the depth of each template based on template amplitude across channels
    templateDepths = list()
    nTemplates = templates.shape[0]
    for templateIndex in np.arange(nTemplates):
        templateEnergyByChannel = np.abs(templates[templateIndex]).sum(0)
        templateDepth = np.average(
            channelDepths,
            weights=templateEnergyByChannel
        )
        templateDepths.append(templateDepth)
    templateDepths = np.array(templateDepths)

    # Assign spike depths based on their corresponding templates
    spikeDepths = templateDepths[templateIndices]

    return templateDepths, spikeDepths

def defineTransformationMatrices(resolution=10):
    """
    Compute the transform that converst CCF voxels to stereotaxic coordinates
    (relative to bregma)
    """

    # Define transformation constants
    bregma = np.array([520, 570.5, 44])
    scale = np.array([-1.031, 0.952, 0.885]) / (1000 / resolution)
    theta = np.radians(5)  # Rotation angle in radians

    # Translation matrix to center at bregma
    T = np.eye(4)
    T[:3, 3] = -np.array(bregma)

    # Scaling matrix to adjust units and reflection
    S = np.diag(np.concatenate([scale, [1]]))

    # Rotation around the y (ML) axis
    R = np.array([
        [np.cos(theta), 0, -np.sin(theta), 0],
        [0,             1,  0,             0],
        [np.sin(theta), 0,  np.cos(theta), 0],
        [0,             0,  0,             1],
    ])

    # Define the forward (F) transformation
    global CCF_TO_STC
    CCF_TO_STC = T @ S @ R
    CCF_TO_STC = R @ S @ T

    # Define the inverse (B) transformation
    global STC_TO_CCF
    R = R.T
    S = np.diag(np.concatenate([1 / scale, [1]]))
    T[:, 3] *= -1
    STC_TO_CCF = T @ S @ R

    return

def transform(points, source='ccf'):
    """
    Apply transformation to coordinates

    Keywords
    --------
    source: str
        Flag which specifies the source space, 'ccf' for the common coordinate
        framework, or 'stc' for stereotaxic coordinates
    """

    #
    if source == 'ccf':
        global CCF_TO_STC
        if CCF_TO_STC is None:
            defineTransformationMatrices()
        tform = CCF_TO_STC
    elif source == 'stc':
        global STC_TO_CCF
        if STC_TO_CCF is None:
            defineTransformationMatrices()
        tform = STC_TO_CCF

    #
    nPoints = points.shape[0]
    homogenousPoints = np.hstack((points, np.ones((nPoints, 1))))

    # Apply the affine transformation matrix
    transformedHomogenousPoints = homogenousPoints @ tform.T

    # Extract the transformed [ML, AP, DV] coordinates
    transformedPoints = transformedHomogenousPoints[:, :3]

    return transformedPoints

def localizeUnits(
    workingDirectory,
    electrodePositionFile=None,
    insertionPoint=None,
    insertionDepth=None,
    insertionAngle=None,
    skullThickness=0.3,
    resolution=10,
    ):
    """
    """

    #
    if type(workingDirectory) == str:
        workingDirectory = pl.Path(workingDirectory)

    #
    if electrodePositionFile is not None:
        electrodePosition = io.loadmat(electrodePositionFile)['probe_positions_ccf'][0][0]
        B = electrodePosition[:, 0]
        A = electrodePosition[:, 1]

    else:
        if insertionPoint is None:
            raise Exception()
        if insertionDepth is None:
            raise Exception()
        if insertionAngle is None:
            raise Exception()
        A = np.array(insertionPoint)
        A[2] += skullThickness
        B = solveForElectrodeEndpoint(
            A,
            insertionDepth,
            insertionAngle
        )
        AB = np.vstack([A, B])
        AB = transform(AB, source='stc')
        A, B = map(np.ravel, np.split(AB, 2, axis=0))

    #
    electrodeVector = A - B
    scalingFactor = electrodeVector / np.linalg.norm(electrodeVector)

    # Load the spike clusters data
    spikeClustersFile = workingDirectory.joinpath('spike_clusters.npy')
    if spikeClustersFile.exists():
        spikeClusters = np.load(spikeClustersFile)

    # Estimate the distance from the tip of the electrode for each unit
    templateDepths, spikeDepths = estimateSpikeDepths(workingDirectory)

    # CCF coordinates for each unit with shape N units x 3 (AP, DV, ML)
    unitCoordinates = list()
    for spikeCluster in np.unique(spikeClusters):
        spikeDepth = np.mean(spikeDepths[spikeClusters == spikeCluster])
        unitCoordinate = B + ((spikeDepth / resolution) * scalingFactor)
        unitCoordinates.append(unitCoordinate)
    unitCoordinates = np.array(unitCoordinates)

    # Lookup the brain area associated with each coordinate
    voxelIndices = np.around(unitCoordinates, 0).astype(int)
    brainStructureIdentities = list()
    global ANNOTATIONS
    if ANNOTATIONS is None:
        _loadAnnotations()
    for (i, j, k) in voxelIndices:
        id = ANNOTATIONS[i, j, k]
        brainStructureIdentities.append(id)
    brainStructureIdentities = np.array(brainStructureIdentities)

    # Convert from CCF voxels to stereotaxic coordinates
    unitCoordinatesTransformed = transform(
        unitCoordinates,
        source='ccf'
    )

    return brainStructureIdentities, unitCoordinates, unitCoordinatesTransformed
