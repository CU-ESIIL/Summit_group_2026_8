import sys
import torch
import torch_geometric
import pandas
import numpy
import sklearn

print("python:", sys.executable)
print("torch:", torch.__version__)
print("torch cuda build:", torch.version.cuda)
print("cuda available:", torch.cuda.is_available())
print("device count:", torch.cuda.device_count())
print("device:", torch.cuda.get_device_name(0) if torch.cuda.is_available() else "none")
print("pyg:", torch_geometric.__version__)
print("pandas:", pandas.__version__)
print("numpy:", numpy.__version__)
print("sklearn:", sklearn.__version__)

if torch.cuda.is_available():
    x = torch.randn(1000, 1000, device="cuda")
    y = x @ x
    print("cuda matmul ok:", y.shape, y.device)
