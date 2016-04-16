#! /bin/sh
# cf. http://pages.uoregon.edu/noeckel/MakeMovie.html

# first convert an image sequence to a movie
ffmpeg -sameq -i %03d.jpg output.mp4

# ... and then convert the movie to a GIF animation
ffmpeg -i output.mp4 -pix_fmt rgb24 -s qcif -loop_output 0 output.gif