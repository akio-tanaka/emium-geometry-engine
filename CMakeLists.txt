cmake_minimum_required(VERSION 3.10)
project(emium-geometry-engine)

# C++の標準を設定
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# libiglのオプション設定
option(LIBIGL_WITH_OPENGL "Use OpenGL" ON)
option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW" ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui" ON)

# サブディレクトリの追加
add_subdirectory(external/libigl)

# gladのインクルード
include(external/libigl/cmake/recipes/external/glad.cmake)

# glfwのインクルード
include(external/libigl/cmake/recipes/external/glfw.cmake)

# インクルードディレクトリの設定
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${LIBIGL_INCLUDE_DIRS}
)
include_directories(${CMAKE_SOURCE_DIR}/external/json/include)

# ソースファイルの収集
file(GLOB SRC_FILES ${CMAKE_SOURCE_DIR}/src/*.cpp)

# 実行ファイルの生成
add_executable(${PROJECT_NAME} ${SRC_FILES})

# ライブラリのリンク
target_link_libraries(${PROJECT_NAME} igl::core glad::glad glfw::glfw)
