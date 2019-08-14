SET(ConstitutiveLaws_src
    ${CMAKE_CURRENT_SOURCE_DIR}/source/BaseMaterial.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/CompressibleNeoHookean.cpp)

add_library( ConstitutiveLaws SHARED ${ConstitutiveLaws_src} )

target_include_directories(
    ConstitutiveLaws
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

DEAL_II_SETUP_TARGET( ConstitutiveLaws )

TARGET_LINK_LIBRARIES( ConstitutiveLaws ${Boost_LIBRARIES} )