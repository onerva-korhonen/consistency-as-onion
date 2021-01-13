# consistency-as-onion
Code for modelling the human brain function with a multi-layer network with flexible nodes. In this network, each layer corresponds to a time window. Nodes on each layer are Regions of Interest (ROIs) optimized in terms of spatial consistency. Intra-layer links are defined as correlation coefficients between ROI time series while inter-layer links indicate the number of shared voxels between ROIs in two layers.

NOTE: I don't update this project on regular basis. For the current state of the flexible-node model, see https://github.com/ercco/multilayer-brains.

*functions.py* contains the functions for creating the flexible-node multilayer as well as for validation of the model. *onion_parameters.py* contains the parameters of the present project. *tests.py* contains a set of tests for the functions; if everything is OK, running this file should return only positive prints.

*atlases* contains possible ROI atlases that could be used as starting points of the ROI optimization (at the moment, there is only the Brainnetome atlas, Fan *et al.* 2016).

*nifti-tools/nifti-functions* contains tools for reading, writing, and manipulating NIFTI files. The NiBabel package (http://nipy.org/nibabel/index.html) plays an important role in these tools.

*scripts* contains a bunch of frontend scripts related to the present project.

*QA.txt* is a list of random mostly programming-related questions that have come to my mind during the project. Unfortunately, I don't promise it to be useful for the larger audience (but don't hesitate to have a look, if you want). *README.md* is the file you are reading right now.

References:

Fan, L., Li, H., Zhuo, J., Zhang, Y., Wang, J., Chen, L., . . .  Jiang, T. (2016). The human Brainnetome atlas:  A new brain atlas based on connectional architecture. *Cerebral Cortex*, 26(8), 3508–3526.


