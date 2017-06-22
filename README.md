# Joint-Spectral-Correspondence
Implementation of Joint Spectral Correspondence for matching the images with disparate appearance arising from factors like dramatic illumination (day vs. night), age (historic vs. new) and rendering style differences.

## Project Description
- Graph Spectral Analysis of image has been used for classic problems like image segmentation, classification, clustering. This project aims to match two images of the same place with dramatic lighting differences using joint spectral analysis of both image.
- This method is based on ![this](Bansal_Joint_Spectral_Correspondence_2013_CVPR_paper.pdf) research paper accepted in CVPR-2013.

- ![src/eigen_function](src/eigen_function) has code to compute joint eigen functions of both images.
- ![src/matching](src/matching) has code to match eigen functions of both images.

- ![dataset](dataset) contains two images for each place and homography between those two images.
- ![results](results) contains precomputed eigen functions in mat file for each image pair of dataset. Matching can be done by directly loading those mat files.
- ![reports](reports) contains following :
	- Theory-basic of eigen spectrum of image : ![image-graph-spectrum-theory.pdf](reports/image-graph-spectrum-theory.pdf)
	- Implementation details of this project : ![implementation.pdf](reports/implementation.pdf)
	- Generated eigen functions for some image pairs of dataset : ![eigen.pdf](reports/eigen.pdf)
	- Results after matching eigen functions of some image pairs : ![matching.pdf](reports/matching.pdf)
