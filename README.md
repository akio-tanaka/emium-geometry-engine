# emium geometry engine

## development

### install dependencies

#### libigl

```sh
git submodule add https://github.com/libigl/libigl.git external/libigl
git submodule update --init --recursive
```

#### nlohmann

```sh
git submodule add https://github.com/nlohmann/json.git external/json
git submodule update --init --recursive
```

#### other library (only for Linux)

```sh
sudo apt update
sudo apt install libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev gdb
```

### build

#### Windows

```sh
mkdir build
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```
- remove `--config Release` when you debug the code.

#### Linux

```sh
mkdir build
cmake -S . -B build
cmake --build build --config Release
```
- remove `--config Release` when you debug the code.

### debug


## usage