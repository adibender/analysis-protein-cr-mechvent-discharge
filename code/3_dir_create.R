#***************************** dir_create *************************************#
#* STEP 3
#* !!! IMPORTANT:  RUN this script only ONCE during the first run
#* !! Treat this carefully if run on  FAT/FAT32 or network-mounted file system
#* This file set up the necessary directories needed for the output
#*
#******************************************************************************#
#*
#***************************** Parameters *************************************#
# Set the name for main output folder used
main_subfolder <- "output"

#******************************************************************************#
# Creates the output directory, to be used as a main directory for the following
dir.create(main_subfolder, showWarnings = FALSE)

# Creation of sub-folders:
dir.create(paste0(main_subfolder,"/data"), showWarnings = FALSE)
dir.create(paste0(main_subfolder,"/models"), showWarnings = FALSE)
dir.create(paste0(main_subfolder,"/results"), showWarnings = FALSE)
dir.create(paste0(main_subfolder,"/results/main"), showWarnings = FALSE)
dir.create(paste0(main_subfolder,"/results/sens"), showWarnings = FALSE)
dir.create(paste0(main_subfolder,"/results/subgroups"), showWarnings = FALSE)
dir.create(paste0(main_subfolder,"/results/figures"), showWarnings = FALSE)
dir.create(paste0(main_subfolder,"/results/figures/supplement"), showWarnings = FALSE)
