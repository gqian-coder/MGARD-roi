add_subdirectory(ParallelHuffman)
list(APPEND MGARD_X_HEADER
    ${CMAKE_CURRENT_SOURCE_DIR}/LZ4.hpp
    )
set(MGARD_X_HEADER ${MGARD_X_HEADER} PARENT_SCOPE)