# TP Particules — Simulation de dynamique moléculaire

Simulation 2D d'impact d'un disque de Lennard-Jones sur une membrane, par méthode des cellules liées et schéma d'intégration de Störmer–Verlet.

## Structure

```
.
├── CMakeLists.txt
├── README.md
├── doc/        Documentation Doxygen + analyse de profiling
├── include/    Déclarations (.hxx)
├── src/        Implémentations (.cxx)
├── test/       Tests unitaires GoogleTest
└── demo/       Programme de démonstration (impact disque/membrane)
```

## Dépendances

- CMake ≥ 3.16
- Compilateur C++17 (GCC, Clang)
- Doxygen et Graphviz (pour `make doc`)
- Paraview (pour visualiser les sorties VTK)

```bash
sudo apt install cmake g++ doxygen graphviz paraview
```

GoogleTest est récupéré automatiquement par CMake (`FetchContent`).

## Compilation

### Mode standard (optimisations max)

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### Mode profiling (gprof)

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Profile ..
make -j$(nproc)
```

## Lancer la simulation

```bash
cd build
./demo/demo
```

Les frames VTK sont écrites dans `build/frames/`. Visualisation avec :

```bash
paraview build/frames/frame_..vtu
```

## Tests

```bash
cd build
ctest --output-on-failure
```

## Documentation

```bash
cd build
make doc
```

La documentation HTML est générée dans `doc/html/`. Ouvrir `doc/html/index.html` dans un navigateur.

## Profiling

Compiler en mode `Profile`, lancer la démo (réduire `tend` dans `demo/demo.cxx` au préalable), puis :

```bash
gprof build/demo/demo build/gmon.out > doc/analyse.txt
```

## Auteurs

Ayoub Touati & Mouad Ikne — Ensimag MMIS 2A