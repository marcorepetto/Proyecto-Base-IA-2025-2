# Pintando con Algoritmos 

¡Hola!  Este repositorio contiene lo necesario para **pintar imágenes en un lienzo** y **convertir instancias en un lienzo**, de manera que ustedes se centren en lo importante: **mejorar las pinceladas** como tal.

Cuentan con los siguientes archivos básicos:

- `stroke.h` — Definiciones de **Canvas**, **ImageGray**, **Stroke** y firmas de utilidades.  
- `stroke.cpp` — Implementación de dibujo, carga/guardado de imágenes y render (pintar). No hace falta que lo entiendas todo.
- `testCall.cpp` — Ejemplo de uso: carga brushes, crea strokes, renderiza (pinta) y guarda un PNG. ¡Compila y ve el ejemplo!

## Conceptos básicos
Un poco de vocabulario en común para entendernos entre todos :).

### 1. Canvas o lienzo
El lienzo es una imagen RGB aplanada de tamaño `width * height * 3` bytes: `[R0, G0, B0, R1, G1, B1, ...]`. Esto lo hace más eficiente para recorrer y guardar.

Por defecto se **inicializa** en blanco, y cuenta con los siguientes helpers:

* `clear(r,g,b)` — pinta todo el lienzo de un color.  
* `setPixel(x,y,r,g,b)` — pinta un píxel (con validación de bordes). Recibe la posición `(x,y)` y el color RGB.  

#### ¿Cómo accedo píxel a píxel?
Cada píxel se guarda en el vector `rgb` de forma lineal, en grupos de 3 valores (R, G, B).  

Para acceder al píxel en `(x,y)` (con `0 ≤ x < canvas.width` y `0 ≤ y < canvas.height`) se usa:

```cpp
int idx = (y * canvas.width + x) * 3; // índice base (R) en el arreglo aplanado
uint8_t R = canvas.rgb[idx + 0];
uint8_t G = canvas.rgb[idx + 1];
uint8_t B = canvas.rgb[idx + 2];

std::cout << "Píxel en (" << x << "," << y << ") tiene el color: "
          << "R=" << (int)R << " G=" << (int)G << " B=" << (int)B << "\n";

```

### 2. Brushes o pinceles
Los **brushes** son imágenes en escala de grises que se usan como **máscaras**:  
- `0` = transparente  
- `255` = opaco  

El color del stroke se aplica modulando por esa máscara.  

* Se cargan con `loadImageGray(filepath, ImageGray&)`.  
* Se guardan en el vector global `gBrushes`.  
* Cada `Stroke` referencia un `type` (índice en `gBrushes`, que va entre 0 y 3).  

### 3. Stroke o pinceladas
Cada pincelada es una estructura que contiene los siguientes atributos:

* `x_rel`, `y_rel` ∈ [0.0,1.0]: centro del stroke en coordenadas relativas.  
* `size_rel` ∈ [0.0,1.0]: tamaño relativo al lado menor del canvas.  
* `rotation_deg` rotation_deg ∈ [0.0 , 360.0): ángulo de rotación en grados (antihorario).
* `type` ∈ [0,3]: índice del brush a usar (`gBrushes[type]`).
* `r,g,b`: color del stroke (cada componente ∈ [0,255]).  

Además, cuenta con un método implementado, que pinta el stroke en el canvas:

* `draw(Canvas& C) const` —  
  - Agarra la **máscara en gris** (brush) que se señala en `type`.  
  - La **escala** al tamaño del stroke (`size_rel`).  
  - La **posiciona** en el lienzo (`x_rel`, `y_rel`).  
  - La **rota** en el lienzo según `rotation_deg`.  

Por cada píxel hace un **muestreo bilineal** de la máscara para suavizar los bordes. Finalmente, mezcla (`alpha-blend`) el color del stroke con el color del fondo, según la intensidad de la máscara. En simples palabras: el stroke pinta el lienzo con su color, pero solo en las zonas que la máscara indica, y con la intensidad que la máscara dicta.

### ¡Importante!
**¿Por qué usar coordenadas relativas?**  Para que los strokes no dependan de un tamaño fijo; cambiar el lienzo a otra resolución mantiene el layout proporcional. _(En caso de que quieran aventurarse con lienzos más grandes)_

## 4. Funciones útiles

* `render(const std::vector<Stroke>&, Canvas&)` —  Recibe un lienzo y un vector de strokes.  
  - Limpia el lienzo a blanco.  
  - Dibuja los strokes **en el orden en que aparecen en el vector**.  

* `savePNG(const Canvas& C, const std::string& filename)` — Guarda un lienzo (`Canvas`) como archivo PNG de 8 bits, y así puedan ver la obra de arte que crean en cada ejecución :).  
  - `C`: el lienzo a guardar.  
  - `filename`: nombre del archivo destino (ej. `"output.png"`).   

* `loadImageRGB_asCanvas(const std::string& filename, Canvas& out)` — Carga una imagen RGB desde el PC (JPG/PNG) y la convierte en un `Canvas`.  
  - `filename`: ruta de la imagen (ej. `"instancias/mona.png"`).  
  - `out`: referencia donde se guarda el canvas resultante.  

Esto será muy útil cuando quieran **comparar** la imagen generada con una **instancia o target**. Así tendrán dos `Canvas` (el generado y el target) y podrán aplicar el MSE para evaluar la calidad de la solución.  

Ejemplo de uso:  
```cpp
Canvas target(1,1); //1px * 1px simplemente para incializar, luego cambia su tamaño automáticamente
if (loadImageRGB_asCanvas("instancias/mona.png", target)) {
    std::cout << "Imagen cargada: " 
              << target.width << "x" << target.height << "\n";

    // Acceder al píxel del centro - ¡Lo que vimos en la sección de lienzo!
    int x = target.width / 2;
    int y = target.height / 2;

    int idx = (y * target.width + x) * 3; // 
    uint8_t R = target.rgb[idx + 0];
    uint8_t G = target.rgb[idx + 1];
    uint8_t B = target.rgb[idx + 2];

    std::cout << "Pixel en (" << x << "," << y << "): "
              << "R=" << (int)R << " G=" << (int)G << " B=" << (int)B << "\n";
}
```

## Para finalizar

Este código base es **tu punto de partida**: puedes modificar lo que estimes necesario según tu enfoque. Sin embargo, no pierdas de vista el objetivo: mejorar las pinceladas.

¿Qué viene después? **Todo**. Este código solamente **lee instancias**, las transforma a un formato comparable con tu imagen generada y convierte tu resultado en algo visible. A partir de aquí debes implementar la parte algorítmica para optimizar el resultado: la creación de la solución inicial, el cálculo del MSE, los movimientos, etc.

### Créditos

- **Brushes**: tomados del repositorio de Anastasia opara’s “Genetic Drawing”. Refererencia [1] del paper de referencia.

- **Instancias**:
  - `mona.png` — *Leonardo da Vinci*, **La Gioconda (Mona Lisa)**. Resolución de 48x64 píxeles.   
  - `dali.png` — *Salvador Dalí*, **La persistencia de la memoria**. Resolución de 64x48 píxeles.   
  - `klimt.png` — *Gustav Klimt*, **El beso**. Resolución de 48x64 píxeles. 
  - `mondriaan.png` — *Piet Mondrian*, **Composición en rojo, azul y amarillo**. Resolución de 48x64 píxeles. 
  - `pollock.png` — *Jackson Pollock*, **Convergence**. Resolución de 64x48 píxeles. 
  - `starrynight.png` — *Vincent van Gogh*, **La noche estrellada**. Resolución de 64x48 píxeles. 
  - `bach.png` — *Elias Gottlob Haussmannm*, **Retrato de Johann Sebastian Bach**. Resolución de 48x64 píxeles. 

- **Librerías externas**:
  - `stb_image.h` y `stb_image_write.h` de Sean Barrett. Licencia public domain / MIT.


