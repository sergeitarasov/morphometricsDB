
#   ____________________________________________________________________________
#   Digitize Landmarks with StereoMorph                                     ####

# To do:
# 1. scale images to W=2000 in GIMP (algorithms do not work with more pixels)
# 2. NOTE min two landmarks should be set for each image
# 3. StereoMorph reference files for landmarks and curves are landmarks_ref.txt and curves_ref.txt
# 4. Images are in data/im
# 5. Ldks are in stereo_output
# 6. Click N as a shortcut for jumping between landmarks

digitizeImage(image.file='data/im', shapes.file='data/stereo_ldk/stereo_output',
              landmarks.ref='data/stereo_ldk/landmarks_ref.txt'
               ) #curves.ref='data/curves_ref.txt'


#   ____________________________________________________________________________
#   Import Landmarks and Save them as RDS                                   ####

rm(ldk)
ldk <- importSM(path='data/stereo_ldk/stereo_output', class = c('Ldk'), panel = TRUE)
str(ldk)
Momocs::get_ldk(ldk)
stereo.ldk <- Momocs::get_ldk(ldk)

# write imported ldks to file
saveRDS(stereo.ldk, 'data/stereo_ldk/stereo.ldk.RDS')
