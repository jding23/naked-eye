# Multi-screen splicing display model

*[HoloLab](http://www.holoddd.com/), Tsinghua University*


<p align="left">
<img src="imgs/principle/see.png">
</p>

<p align="left"> <strong>Figure 1</strong>.The multi-screen splicing model uses a perspective view, where the position of the viewpoint is very important. The relative position of the viewpoint to different screens is different, and parameters need to be transformed during the mapping process to meet the requirements of the stereoscopic effect. </p>

This model primarily achieves a 3D presentation effect of point cloud maps on the splicing screen. It mainly uses a layer-by-layer loading method for mapping point cloud projections onto the screen, segmenting the point cloud, and constructing perspective mapping relationships for each segment.

## Reference
Matlab 2019b or newer. Older visions may be sufficient but have not been tested.

## Quick Start
- **Splicing screen effect demonstration.** [`demo`](https://github.com/jding23/naked-eye/edit/main/coding/demo) includes two types of demonstrations: dual-screen splicing and three-screen splicing.
- **Parameter Settings:** The screen parameters used refer to the existing products of Huawei    [`two_screens_display.mat`](https://github.com/jding23/naked-eye/edit/main/coding/config/two_screens_display.mat) and Samsung [`three_screens_display.mat`](https://github.com/jding23/naked-eye/edit/main/coding/config/three_screens_display.mat).
- **Point cloud:** The Stanford bunny point cloud model is selected as the experimental object[`bun_zipper.ply`](https://github.com/jding23/naked-eye/edit/main/coding/point_cloud/bun_zipper.ply).It has no color.
- **Background:** Remove the screen face from the six-sided sky box.

## Theories and References
Please refer to the article[``]