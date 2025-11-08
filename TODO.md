Objetivo general: generar una configuración inicial de pinceladas (posición, orientación, tamaño) que capture la **estructura geométrica** de la imagen objetivo (bordes, texturas, contornos dominantes) sin preocuparse aún del color.
Esa configuración estructural servirá como base para la optimización de color posterior por mínimos cuadrados.

En tu implementación actual, el `initialGuess()` evalúa cada pincelada candidata midiendo si cubre zonas con bordes (`targetImageBorders.gray`), penalizando zonas ya cubiertas.
Eso equivale a un criterio binario: “¿cuánto borde cubrí o no cubrí?”.
El problema es que este criterio **no distingue entre tipos de detalle** (bordes gruesos vs. finos, textura vs. superficie plana), ni considera la **orientación local del borde**.

Queremos reemplazar esa heurística limitada por una **función de utilidad estructural** que guíe las pinceladas hacia zonas de alta energía espacial (bordes y texturas relevantes) y evite redundancia.

---

### Qué estamos intentando lograr

1. **Cubrir información estructural significativa**:
   Las primeras pinceladas deben ubicarse donde la imagen tiene mayor información visual: bordes fuertes, patrones finos, líneas de contorno.

2. **Evitar redundancia**:
   Si varias pinceladas cubren la misma región, el beneficio marginal cae; debemos penalizar solapamientos innecesarios.

3. **Proporcionar base para el color**:
   Si la estructura está bien distribuida espacialmente, la etapa de color podrá interpolar mejor —sin necesidad de mover pinceladas después.

4. **Reducir el error estructural global (E(x,y))**:
   La función greedy debe buscar maximizar la reducción de la energía estructural total de la imagen con cada nueva pincelada.

---

### Qué modificaciones introducimos en `initialGuess()`

1. **Construcción del mapa estructural (E(x,y))**
   Antes del bucle de pinceladas, calculamos un mapa que combine:

   * **Magnitud del gradiente** (bordes):
     (G(x,y) = |\nabla I(x,y)|)
   * **Textura local** (variación local o energía de alta frecuencia):
     (T(x,y) = \text{Var}_{\text{local}}(I))

   Luego los combinamos:
   [
   E(x,y) = \alpha , G(x,y) + \beta , T(x,y)
   ]
   Este mapa representa cuánta información estructural hay en cada píxel.

2. **Nuevo criterio greedy para seleccionar cada trazo**
   En lugar de usar la métrica de cobertura binaria (`covered_pixels` y `uncovered_pixels`), medimos cuánto “ganamos” en términos de reducción de energía estructural:
   [
   \text{gain} = \sum_{p \in \text{pincelada}} a_p , E(p)
   ]
   y penalizamos el solapamiento con trazos previos:
   [
   \text{penalty} = \sum_{p \in \text{pincelada}} a_p , C(p)
   ]
   donde (C(p)) es el mapa acumulado de cobertura.
   El puntaje final del trazo:
   [
   \text{score} = \text{gain} - \lambda , \text{penalty}
   ]

3. **Actualización del mapa (E)**
   Cuando se selecciona una pincelada, la energía estructural en los píxeles que cubre se reduce:
   [
   E(p) \leftarrow E(p) \cdot (1 - a_p)
   ]
   Esto evita que los siguientes trazos sigan eligiendo las mismas zonas.

4. **Resultado esperado**

   * Las primeras pinceladas se ubican sobre los bordes más fuertes.
   * Las siguientes se distribuyen en texturas y detalles residuales.
   * Las zonas planas quedan para el ajuste de color posterior.
   * Se obtiene una base geométrica más informativa y estable.

---

### Resumen conceptual

| Componente               | Estado actual                             | Modificación                                              |
| ------------------------ | ----------------------------------------- | --------------------------------------------------------- |
| Métrica de evaluación    | Basada en cobertura binaria de bordes     | Basada en energía estructural continua (bordes + textura) |
| Penalización             | Cobertura simple                          | Penalización proporcional a solapamiento                  |
| Actualización del estado | Incrementa mapa de bordes                 | Reduce energía estructural residual                       |
| Efecto final             | Pinceladas concentradas en bordes gruesos | Pinceladas adaptadas a estructura global y textura local  |

---

En resumen:
**lo que intentamos hacer** es redefinir la lógica de “qué significa una buena pincelada inicial”.
Dejar de medir solo “cuánto borde cubro” y pasar a “cuánta información estructural útil capturo y cuán poco redundante soy”.
Esa información estructural luego define la geometría base que la optimización de color podrá usar eficientemente.
