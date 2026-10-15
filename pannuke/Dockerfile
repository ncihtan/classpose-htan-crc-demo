FROM pytorch/pytorch:2.5.1-cuda12.1-cudnn9-runtime

RUN apt-get update && apt-get install -y --no-install-recommends git && rm -rf /var/lib/apt/lists/*

# Install classpose and training deps (bypassing Python>=3.13 version pin;
# classpose code itself is compatible with Python 3.11)
RUN pip install --no-cache-dir --ignore-requires-python \
    "classpose @ git+https://github.com/sohmandal/classpose@main" \
    datasets scikit-learn

# classpose upgrades torch to cu128; reinstall cu121 to match the system CUDA
# (pytorch/pytorch base has system CUDA 12.1; requires driver >=525, well supported on A100)
RUN pip install --no-cache-dir --force-reinstall \
    "torch==2.5.1+cu121" "torchvision==0.20.1+cu121" \
    --index-url https://download.pytorch.org/whl/cu121

# Patch missing transforms subpackage
RUN python - <<'EOF'
import os, urllib.request, importlib.util
spec = importlib.util.find_spec("classpose")
pkg_dir = os.path.dirname(spec.origin)
transforms_dir = os.path.join(pkg_dir, "transforms")
os.makedirs(transforms_dir, exist_ok=True)
COMMIT = "9d18066bb4e5b7d95e00a7d913b42663e228dc9b"
BASE = f"https://raw.githubusercontent.com/sohmandal/classpose/{COMMIT}/src/classpose/transforms"
for fname in ["__init__.py","augmentation_configs.py","he_staining.py","hed.py","image_quality.py","stardist_augmentation.py","transforms.py"]:
    urllib.request.urlretrieve(f"{BASE}/{fname}", os.path.join(transforms_dir, fname))
    print(f"patched: {fname}")
EOF

WORKDIR /app
COPY run_training.py .
COPY trainer/ trainer/

ENTRYPOINT ["python", "-m", "trainer.task"]
