
# The text files here are used to use image making
# They need to be found via FW_SEARCH_PATH

# install constraint files
file(GLOB constraints *.txt)
install_fw(LIST ${constraints})
