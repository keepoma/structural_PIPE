import numpy as np
import matplotlib.pyplot as plt

# Loading matrix
# The matrix is 379x379: 180 cortical regions on each hemisphere + 19 subcortical
sc = np.genfromtxt('/media/nas/nikita/sample_598/raw/9_atlas/hcpmmp1.csv', delimiter=',')

# Removing maximum element value to fix color coding
upper_bound = np.percentile(sc, 99)  # 99th percentile to avoid extreme outliers
sc_clipped = np.clip(sc, None, upper_bound)

# Display the SC matrix
plt.figure(figsize=(8, 6))
plt.imshow(sc_clipped, cmap='viridis', aspect='auto')
plt.colorbar(label='SC Value')
plt.title('Structural Connectivity Matrix')
plt.show()
