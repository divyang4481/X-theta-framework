param(
  [string]$EnvName = "xtheta"
)

conda env create -f env\environment.yml
conda activate $EnvName

python -m pip install --upgrade pip
python -m pip install --index-url https://download.pytorch.org/whl/cu124 `
  torch==2.6.0+cu124 torchvision==0.21.0+cu124 torchaudio==2.6.0+cu124

python - << 'PYCODE'
import torch
print("Torch:", torch.__version__)
print("CUDA available:", torch.cuda.is_available())
print("CUDA version:", torch.version.cuda)
PYCODE
