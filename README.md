# Description
This is a Python package for estimating the location of cells recorded with Neuropixels within the mouse brain

# Installation
Clone the repository:<br />
`git clone https://github.com/jbhunt/nptracer.git`<br />
<br />
Navigate to the root directory of the project:<br />
`cd ./nptracer`<br />
<br />
Execute the setup script with pip:<br />
`pip install .`

# Basic usage
## Estimating unit location with the Neuropixels trajectory explorer
```Python
import nptracer as npt
kilosortOutputFolder = 'path/to/kilosort/output'
trajectoryExplorerFile = 'path/to/trajectory/explorere/file'
labels, points, transformed = npt.localizeUnits(
    kilosortOutputFolder=kilosortOutputFolder
    trajectoryExplorerFile=trajectoryExplorerFile
)
```
`labels` is an array of the brain structures associated with each unit. `points`
is an array (N units x 3) with the coordinates of each unit in voxels of the
Allen Common Coordinate Framework. `transformed` is an array with the same shape
as `points` but indicates the location of each unit in stereotaxic coordinates
in reference to bregma.<br />
## Estimating unit location with sterotaxic coordinates
```Python
import nptracer as npt
kilosortOutputFolder = 'path/to/kilosort/output'
insertionPoint = np.array([-3.9, 2.5, 0.4]) # Insertion point in sterotaxic coordinates (AP, ML, DV) in mm
insertionDepth = 3.6 # Depth of insertion along the trajectory of the insertion (in mm)
insertionAngle = 6 # Angle of insertion (in degrees)
skullThickness = 0.3 # Assumed thickness of the skull (in mm)
labels, points, transformed = npt.localizeUnits(
    kilosortOutputFolder=kilosortOutputFolder,
    insertionPoint=insertionPoint,
    insertionDepth=insertionDepth,
    insertionAngle=insertionAngle,
    skullThickness=skullThickness
)
```
## Utility functions
```Python
import nptracer as npt
kilosortOutputFolder = 'path/to/kilosort/output'
templateDepths, spikeDepths = npt.estimateSpikeDepths(
    kilosortOutputFolder=kilosortOutputFolder
)
```
This function computes the distance of each spike from the tip of the Neurpixels
probe. Please note that this function assumes all active channels are clustered
at the tip of the electrode.