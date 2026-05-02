::: mermaid
flowchart TD
    A[Inicio] --> B[Importar librerías y definir constantes]
    B --> C[Definir lista de índices n]
    C --> D[Para cada n en la lista]
    D --> E["integrar_lane_emden(n)"]
    E --> F["Inicializar con Taylor"]
    F --> G{"¿theta nueva < 0?"}
    G -- No --> H[Paso RK4, almacenar punto]
    H --> G
    G -- Sí --> I["find_surface: localizar ξ₁ con spline + Brent"]
    I --> J["compute_parameters: calcular ρc, etc."]
    J --> K[Guardar perfil y parámetros]
    K --> L{"¿Último n?"}
    L -- No --> D
    L -- Sí --> M[Crear tabla de resultados]
    M --> N["Graficar θ(ξ) para todos los n"]
    N --> O{"¿Existe solar_model.txt?"}
    O -- Sí --> P[Cargar datos, interpolar y comparar densidades]
    O -- No --> Q[Saltar comparación solar]
    P --> R[Mostrar gráficos y tabla]
    Q --> S
    R --> S[Fin]
:::