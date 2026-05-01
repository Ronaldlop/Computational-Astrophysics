\`\`\`mermaid
flowchart TD
    A[Inicio] --> B[Importar librerías y definir constantes]
    B --> C[Definir n_values]
    C --> D[Para cada n en n_values]
    D --> E[integrar_lane_emden(n)]
    E --> F[Initializar con Taylor]
    F --> G[¿y1_new < 0?]
    G -- No --> H[Paso RK4, almacenar punto]
    H --> G
    G -- Sí --> I[find_surface: localizar ξ1 con spline + Brent]
    I --> J[compute_parameters: obtener ρc, etc.]
    J --> K[Guardar perfil y parámetros]
    K --> L[¿Último n?]
    L -- No --> D
    L -- Sí --> M[Crear tabla de resultados]
    M --> N[Graficar θ(ξ) para todos los n]
    N --> O[¿Existe solar_model.txt?]
    O -- Sí --> P[Cargar datos, interpolar y comparar densidades]
    O -- No --> Q[Fin sin comparación solar]
    P --> Q
    Q --> R[Mostrar gráficos y tabla]
    R --> S[Fin]
\`\`\`