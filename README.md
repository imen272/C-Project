# Documentation du Projet C++

##  Installation
```bash
# Linux
sudo apt install git cmake g++

# macOS
brew install git cmake
```

##  Lancer le projet
```bash
git clone https://github.com/imen272/C-Project.git
cd C-Project
mkdir -p build && cd build
cmake .. && make
```

## Tests
```bash
ctest --output-on-failure
```
## Lancer Executable
```bash
./demo/demo_univers

```
## 📊 Visualisation (ParaView)
```bash
paraview output/particules_t_*.vtu
```

## Aide
| Commande | Action |
|----------|--------|
| `make clean` | Nettoie la compilation |
| `rm -rf build` | Réinitialise tout |
