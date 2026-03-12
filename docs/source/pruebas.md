# tests - Estrategia de validacion

Este directorio contiene pruebas unitarias y de consistencia para los componentes criticos del modelo.

## Archivos de prueba

- `test_utils.py`: robustez numerica (`_safe_exp`, `_finite_or_zero`).
- `test_lattice.py`: geometria SOS, contorno periodico y reglas de migracion.
- `test_bkl.py`: consistencia de logica BKL, clasificacion y balance de masa.
- `tests_suite.py`: suite integrada de pruebas.

## Cobertura funcional esperada

- Estabilidad numerica en tasas y exponenciales.
- Integridad geometrica de la red.
- Conservacion y consistencia de estados en eventos estocasticos.

## Criterios minimos antes de publicar documentacion

- Sin fallos en pruebas unitarias.
- Sin errores de import en modulos de `src/`.
- Sin warnings criticos al construir Sphinx en modo estricto.

## Notas para Sphinx

Este README puede enlazarse desde la seccion de uso para dejar explicita la politica de validacion tecnica.
