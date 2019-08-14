SET(FE_src
    ${CMAKE_CURRENT_SOURCE_DIR}/source/QPData.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/SimulationParameters.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/RHS_Functions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/SolidSerial.cpp)

add_library( FE SHARED ${FE_src} )

target_include_directories(
    FE
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

DEAL_II_SETUP_TARGET( FE )

TARGET_LINK_LIBRARIES( FE ${Boost_LIBRARIES} )
TARGET_LINK_LIBRARIES( FE ConstitutiveLaws )