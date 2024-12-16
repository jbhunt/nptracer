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

# TODO: Code this function
def downloadCommonCoordinateFrameworkData():
    """
    """

    return

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
    The insertionPoint argument must be a tuple or array like object that
    specifies the stereotaxic coordinates of the insertion point (AP, ML, DV).

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

def _computeTransformMatrix(
    inverse=False
    ):
    """
    Compute the transform that converst CCF voxels to stereotaxic coordinates
    (relative to bregma)
    """

    # Define transformation constants
    bregma = [570.5, 520, 44]  # [ML, AP, DV]
    scale = np.array([0.952, -1.031, 0.885]) / 100  # Scaling factors for [ML, AP, DV]
    rotation = np.radians(5)  # Rotation angle in radians

    # Translation matrix to center at bregma
    T = np.eye(4)
    T[:3, 3] = -np.array(bregma)

    # Scaling matrix to adjust units and reflection
    S = np.diag(np.concatenate([scale, [1]]))

    # Rotation matrix for nose-up tilt
    R = np.array([
        [1, 0, 0, 0],
        [0, np.cos(rotation), -np.sin(rotation), 0],
        [0, np.sin(rotation), np.cos(rotation), 0],
        [0, 0, 0, 1]
    ])

    #
    if inverse:
        R = R.T
        S = np.diag(np.concatenate([1 / scale, [1]]))
        T[:, 3] *= -1
        F = T @ S @ R

    else:
        F = R @ S @ T

    return F

def applyTransform(points, source='ccf'):
    """
    Apply transformation to coordinates

    Keywords
    --------
    source: str
        Flag which specifies the source space, ccf for the common coordinate
        framework, or st for stereotaxic
    """

    # Determine the direction of the transformation
    if source == 'ccf':
        inverse = False
    elif source == 'st':
        inverse = True
    T = _computeTransformMatrix(inverse)

    #
    nPoints = points.shape[0]
    homogenousPoints = np.hstack((points, np.ones((nPoints, 1))))

    # Apply the affine transformation matrix
    transformedHomogenousPoints = homogenousPoints @ T.T

    # Extract the transformed [ML, AP, DV] coordinates
    transformedPoints = transformedHomogenousPoints[:, :3]

    return transformedPoints

def estimateLocationForAllUnits(
    workingDirectory,
    electrodePositionFile,
    resolution=10,
    ):
    """
    Estimate the CCF coordinates for all units in a recording based on electrode
    position and channel mapping
    """

    #
    if type(workingDirectory) == str:
        workingDirectory = pl.Path(workingDirectory)

    # Estimate the distance from the tip of the electrode for each unit
    templateDepths, spikeDepths = estimateSpikeDepths(workingDirectory)

    # Coordinates from the Allen CCF with the shape 3 (AP, DV, ML) x 2 (electrode tip, electrode insertion)
    electrodePosition = io.loadmat(electrodePositionFile)['probe_positions_ccf'][0][0]
    tipPosition = electrodePosition[:, 0]
    entryPosition = electrodePosition[:, 1]
    electrodeVector = entryPosition - tipPosition
    scalingFactor = electrodeVector / np.linalg.norm(electrodeVector)

    # Load the spike clusters data
    spikeClustersFile = workingDirectory.joinpath('spike_clusters.npy')
    if spikeClustersFile.exists():
        spikeClusters = np.load(spikeClustersFile)

    # CCF coordinates for each unit with shape N units x 3 (AP, DV, ML)
    unitCoordinates = list()
    for spikeCluster in np.unique(spikeClusters):
        spikeDepth = np.mean(spikeDepths[spikeClusters == spikeCluster])
        unitCoordinate = tipPosition + ((spikeDepth / resolution) * scalingFactor)
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

    # Re-order axes from AP, DV, ML to ML, AP, DV (necessary for transform)
    unitCoordinates = np.vstack([
        unitCoordinates[:, 2],
        unitCoordinates[:, 0],
        unitCoordinates[:, 1]
    ]).T

    # Convert from CCF voxels to stereotaxic coordinates
    unitCoordinatesTransformed = applyTransform(
        unitCoordinates,
        source='ccf'
    )

    return brainStructureIdentities, unitCoordinates, unitCoordinatesTransformed
