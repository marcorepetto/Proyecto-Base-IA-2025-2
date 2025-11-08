# --- Directorios de Salida ---
# Directorio para los archivos objeto (.o)
OBJ_DIR = obj
# Directorio para el ejecutable final
BIN_DIR = bin

# --- Variables del Compilador ---
CXX = g++
# Opciones de compilación (Incluir headers, advertencias, depuración)
CXXFLAGS = -I ./include/ -Wall -g

# --- Variables de Archivos ---
# El nombre del ejecutable final, ahora DENTRO de bin/
TARGET = $(BIN_DIR)/programa

# La lista de archivos .o, ahora DENTRO de obj/
OBJS = $(OBJ_DIR)/main.o $(OBJ_DIR)/stroke.o $(OBJ_DIR)/model.o

n ?= 125
k ?= 100
pixel_threshold ?= 0.1
p ?= 0.5

# --- Reglas ---

# Regla principal (lo que se ejecuta si solo escribes 'make')
all: $(TARGET)

# Un 'alias' para que 'make compile' también funcione
compile: all

# Regla de enlazado: Cómo crear el ejecutable en 'bin/'
# Depende de todos los archivos .o que están en 'obj/'
$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(OBJS) -o $(TARGET)

# --- Reglas de Compilación (Objeto) ---
# Cada regla ahora debe:
#  1. Asegurar que la carpeta 'obj/' exista
#  2. Usar '-o' para guardar el .o en la carpeta 'obj/'

$(OBJ_DIR)/main.o: main.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c main.cpp -o $(OBJ_DIR)/main.o

$(OBJ_DIR)/stroke.o: stroke.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c stroke.cpp -o $(OBJ_DIR)/stroke.o

$(OBJ_DIR)/model.o: model.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c model.cpp -o $(OBJ_DIR)/model.o

# --- Regla de Limpieza ---
# Borra las carpetas 'obj' y 'bin' enteras
clean:
	@echo "Limpiando carpetas obj/ y bin/..."
	@rm -rf $(OBJ_DIR) $(BIN_DIR)

# --- Reglas Especiales ---
run: $(TARGET)
	./$(TARGET) bach.png $(n) $(k) $(pixel_threshold) $(p)
	./$(TARGET) dali.png $(n) $(k) $(pixel_threshold) $(p)
	./$(TARGET) klimt.png $(n) $(k) $(pixel_threshold) $(p)
	./$(TARGET) mona.png $(n) $(k) $(pixel_threshold) $(p)
	./$(TARGET) mondriaan.png $(n) $(k) $(pixel_threshold) $(p)
	./$(TARGET) pollock.png $(n) $(k) $(pixel_threshold) $(p)
	./$(TARGET) starrynight.png $(n) $(k) $(pixel_threshold) $(p)

# Le dice a 'make' que 'all', 'compile' y 'clean' no son
# nombres de archivos reales, para evitar conflictos.
.PHONY: all compile clean run