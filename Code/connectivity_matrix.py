import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('hcpmmp1 Kopie.csv', header=None)

data = df.to_numpy(dtype=float)
upper_bound = np.percentile(data, 99)
data_clipped = np.clip(data, None, upper_bound)

plt.imshow(data_clipped, cmap='plasma')
plt.colorbar(label='SC Value')
plt.title('Structural Connectivity Matrix')
plt.show()
