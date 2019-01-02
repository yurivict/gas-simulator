rm -f video.avi
ffmpeg -start_number 1 -i image-t=0.%05d.png -vcodec mpeg4 -vb 40M video.avi
