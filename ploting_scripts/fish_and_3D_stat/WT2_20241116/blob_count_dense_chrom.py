import cv2
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from tifffile import imread
from scipy.spatial import distance

# Load the TIFF image
# image = imread('WT1-SDK1_CNTNA2.tif')
image = imread('WT2_CNTNA2_RHEB.tif')

swapped_image = image.copy()
swapped_image[:, :, 0], swapped_image[:, :, 2] = image[:, :, 2], image[:, :, 0]
blue_channel = swapped_image[:, :, 0]
cv2.imwrite('01.raw_blue_channel.png', blue_channel)

##############

# Plot the pixel intensity distribution of the raw blue channel
plt.figure(figsize=(10, 6))
plt.hist(blue_channel.ravel(), bins=256, color='blue', alpha=0.7, density=True)
plt.title('Pixel Intensity Distribution of Raw Blue Channel')
plt.xlabel('Pixel Intensity')
plt.ylabel('Density')
plt.xlim(0, 255)
plt.ylim(0, 0.05)
# Save the histogram plot to a PNG file
plt.savefig('02.raw_blue_channel_intensity_distribution.svg')

# Show the plot (optional)
plt.show()
##########################################
# Extract blue channel
blue_channel[blue_channel <= 50] = 0

cv2.imwrite('03.filtered_blue_channel.png', blue_channel)

# Apply GaussianBlur to reduce noise and improve contour detection
blurred = cv2.GaussianBlur(blue_channel, (5, 5), 0)
cv2.imwrite('04.blurred.png', blurred)

# Perform edge detection
edges = cv2.Canny(blue_channel, 30, 30)

# Save the initial edge-detected image
cv2.imwrite('05.edges_initial.png', edges)

# Apply dilation to fill gaps in edges
kernel = np.ones((2, 2), np.uint8)
dilated = cv2.dilate(edges, kernel, iterations=1)

# Save the dilated edge-detected image
cv2.imwrite('06.edges_dilated.png', dilated)

# Apply closing to fill gaps around nearby edges
closed = cv2.morphologyEx(dilated, cv2.MORPH_CLOSE, kernel)

# Save the processed edge-detected image
cv2.imwrite('07.edges_processed.png', closed)

# Find contours in the processed edge-detected image
contours, _ = cv2.findContours(closed, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

# Sort contours by area in descending order and keep the four largest
contours = sorted(contours, key=cv2.contourArea, reverse=True)[:40]

# Draw the contours on the original image and label them
output = swapped_image.copy()
for i, contour in enumerate(contours):
    cv2.drawContours(output, [contour], -1, (0, 255, 0), 2)
    # Calculate the center of the contour for labeling
    M = cv2.moments(contour)
    if M["m00"] != 0:
        cX = int(M["m10"] / M["m00"])
        cY = int(M["m01"] / M["m00"])
    else:
        cX, cY = 0, 0
    cv2.putText(output, f"#{i + 1}", (cX, cY), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 0, 255), 2)

# Save the output image with contours highlighted and labeled
cv2.imwrite('08.contours_detected.png', output)


##########################################
# Convert to HSV color space for better color segmentation
hsv_image = cv2.cvtColor(swapped_image, cv2.COLOR_BGR2HSV)

# Define color ranges for light blue, red, and green blobs
red_lower = np.array([0, 50, 50])
red_upper = np.array([10, 255, 255])
# green_lower = np.array([40, 50, 50])
# green_upper = np.array([80, 255, 255])
# cyan_lower = np.array([80, 50, 50])
# cyan_upper = np.array([100, 255, 255])
purple_lower = np.array([130, 50, 50])
purple_upper = np.array([160, 255, 255])

# Create masks for each color
mask_red = cv2.inRange(hsv_image, red_lower, red_upper)
# mask_green = cv2.inRange(hsv_image, green_lower, green_upper)
# mask_cyan = cv2.inRange(hsv_image, cyan_lower, cyan_upper)
mask_purple = cv2.inRange(hsv_image, purple_lower, purple_upper)

# Save individual masks to PNG files for visual inspection
cv2.imwrite('09.mask_red.png', mask_red)
cv2.imwrite('10.mask_purple.png', mask_purple)


##########################################
# Find contours for red and green blobs
contours_red, _ = cv2.findContours(mask_red, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
# contours_green, _ = cv2.findContours(mask_green, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
# contours_cyan, _ = cv2.findContours(mask_cyan, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
contours_purple, _ = cv2.findContours(mask_purple, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

# Function to check if a point is inside any of the reference contours
def is_point_inside_contours(point, reference_contours):
    for ref_contour in reference_contours:
        if cv2.pointPolygonTest(ref_contour, point, False) >= 0:
            return True
    return False

# Filter blobs to only include those inside the reference contours
def filter_contours_by_reference(contours, reference_contours):
    filtered_contours = []
    for contour in contours:
        M = cv2.moments(contour)
        if M["m00"] != 0:
            cX = int(M["m10"] / M["m00"])
            cY = int(M["m01"] / M["m00"])
            if is_point_inside_contours((cX, cY), reference_contours):
                filtered_contours.append(contour)
    return filtered_contours

filtered_contours_red = filter_contours_by_reference(contours_red, contours)
# filtered_contours_green = filter_contours_by_reference(contours_green, contours)
# filtered_contours_cyan = filter_contours_by_reference(contours_cyan, contours)
filtered_contours_purple = filter_contours_by_reference(contours_purple, contours)

##########################################
# Draw red and green blobs on the original image
image_with_blobs = swapped_image.copy()
# Label red blobs
for i, contour in enumerate(filtered_contours_red):
    cv2.drawContours(image_with_blobs, [contour], -1, (0, 0, 255), 1)  # Red blobs in red color
    M = cv2.moments(contour)
    if M["m00"] != 0:
        cX = int(M["m10"] / M["m00"])
        cY = int(M["m01"] / M["m00"])
        cv2.putText(image_with_blobs, f"R{i+1}", (cX, cY), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 0, 255), 1)

# # Label green blobs
# for i, contour in enumerate(filtered_contours_green):
#     cv2.drawContours(image_with_blobs, [contour], -1, (0, 255, 0), 1)  # Green blobs in green color
#     M = cv2.moments(contour)
#     if M["m00"] != 0:
#         cX = int(M["m10"] / M["m00"])
#         cY = int(M["m01"] / M["m00"])
#         cv2.putText(image_with_blobs, f"G{i + 1}", (cX, cY), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 255, 0), 1)
#
# # Label cyan blobs
# for i, contour in enumerate(filtered_contours_cyan):
#     cv2.drawContours(image_with_blobs, [contour], -1, (255, 255,0), 1)  # cyan blobs in cyan color
#     M = cv2.moments(contour)
#     if M["m00"] != 0:
#         cX = int(M["m10"] / M["m00"])
#         cY = int(M["m01"] / M["m00"])
#         cv2.putText(image_with_blobs, f"C{i + 1}", (cX, cY), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 255,0), 1)

# Label purple blobs
for i, contour in enumerate(filtered_contours_purple):
    cv2.drawContours(image_with_blobs, [contour], -1, (255, 0, 255), 1)  # purple blobs in purple color
    M = cv2.moments(contour)
    if M["m00"] != 0:
        cX = int(M["m10"] / M["m00"])
        cY = int(M["m01"] / M["m00"])
        cv2.putText(image_with_blobs, f"P{i + 1}", (cX, cY), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 0, 255), 1)

# Save the image with red and green blobs in JPG format
cv2.imwrite('11.image_with_blobs.png', image_with_blobs)

# Display the image with red and green blobs (optional)
plt.figure(figsize=(10, 10))
plt.imshow(cv2.cvtColor(image_with_blobs, cv2.COLOR_BGR2RGB))
plt.title('Image with Blobs')
plt.axis('off')
plt.show()
plt.savefig('12.image_with_blobs.svg')


##########################################
# Function to calculate the distance from each blob to its closest contour line
def calculate_distances(blob_contours, reference_contours):
    distances = []
    for blob_contour in blob_contours:
        for point in blob_contour:
            point = point[0]
            min_distance = float('inf')
            for ref_contour in reference_contours:
                for ref_point in ref_contour:
                    ref_point = ref_point[0]
                    dist = distance.euclidean(point, ref_point)
                    if dist < min_distance:
                        min_distance = dist
            distances.append(min_distance)
    return distances

# Calculate distances for red and green blobs to the blue channel contours
distances_red = calculate_distances(filtered_contours_red, contours)
# distances_green = calculate_distances(filtered_contours_green, contours)
distances_purple = calculate_distances(filtered_contours_purple, contours)
# distances_cyan = calculate_distances(filtered_contours_cyan, contours)

# Create a figure and a set of subplots
fig, axs = plt.subplots(2, 2, figsize=(10, 10))

# Plot the distribution of distances for red blobs
axs[0, 0].hist(distances_red, bins=50, color='red', alpha=0.7, density=True)
axs[0, 0].set_title('Distance of CNTNAP2 to Closest nuclear membrane')
axs[1, 0].set_xlabel('Distance')
axs[1, 0].set_ylabel('Density')

# # Plot the distribution of distances for green blobs
# axs[0, 1].hist(distances_green, bins=50, color='green', alpha=0.7, density=True)
# axs[0, 1].set_title('Distance of SDK1 to Closest nuclear membrane')
# axs[1, 0].set_xlabel('Distance')
# axs[1, 0].set_ylabel('Density')
#
# # Plot the distribution of distances for purple blobs
axs[1, 1].hist(distances_purple, bins=50, color='purple', alpha=0.7, density=True)
axs[1, 1].set_title('Distance of RHBE to Closest nuclear membrane')
axs[1, 0].set_xlabel('Distance')
axs[1, 0].set_ylabel('Density')

# Plot the distribution of distances for cyan blobs
# axs[1, 0].hist(distances_cyan, bins=50, color='cyan', alpha=0.7, density=True)
# axs[1, 0].set_title('Distance of CALN1 to Closest nuclear membrane')
# axs[1, 0].set_xlabel('Distance')
# axs[1, 0].set_ylabel('Density')

plt.savefig('13.distance_distribution_blobs.pdf')
plt.show()
