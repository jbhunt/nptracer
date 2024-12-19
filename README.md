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
## Using the trajectory estimated with the Neuropixels trajectory explorer
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
is an array (N units x 3) with the coordinates of each unit in voxels fo the
Allen Common Coordinate Framework. `transformed` is an array with the same shape
ad `points` but indicates the location of each unit in stereotaxic coordinates
in reference to bregma.