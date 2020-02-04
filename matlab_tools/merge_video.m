%mergeVideo
cd /project/3015069.01/model/figures/videos


video_path = '/project/3015069.01/model/figures/videos/merge_videos.avi'

writerObj = VideoWriter(video_path);
writerObj.FrameRate = 10;

vid = VideoReader('axon_packing_video_esmrmb.avi');

while hasFrame(vid)
    F{1} = readFrame(vid);
    keyboard;
end

vid = VideoReader('zoomed_on_model_esmrmb.avi');

while hasFrame(vid)
    F(end+1) = readFrame(vid);
    keyboard;
end

