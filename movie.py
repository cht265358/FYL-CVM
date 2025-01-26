from moviepy.editor import ImageSequenceClip
import os
import imageio

# Directory containing your PNG files
image_folder = 'phasediagram_elastic'

# Get a sorted list of all image files in the folder
images = sorted([os.path.join(image_folder, img) for img in os.listdir(image_folder) if img.endswith(".png")])

# Create a video from the images
output_video = 'output_video.mp4'
with imageio.get_writer(output_video, fps=1) as writer:  # Adjust fps (frames per second) as needed
    for image in images:
        img = imageio.imread(image)  # Read each image
        writer.append_data(img)  # Add the image to the video

print(f"Video saved as {output_video}")