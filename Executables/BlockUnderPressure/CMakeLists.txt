add_executable( BlockUnderPressure BlockUnderPressure.cc )

target_link_libraries( BlockUnderPressure FE )

DEAL_II_SETUP_TARGET(BlockUnderPressure)