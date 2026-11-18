import os
import pickle
import numpy as np

try:
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.linear_model import LogisticRegression
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

try:
    import torch
    import torch.nn as nn
    import torch.optim as optim
    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False

if HAS_TORCH:
    class DeepHistovarNet(nn.Module):
        def __init__(self, input_size):
            super(DeepHistovarNet, self).__init__()
            self.fc1 = nn.Linear(input_size, 16)
            self.relu1 = nn.ReLU()
            self.fc2 = nn.Linear(16, 8)
            self.relu2 = nn.ReLU()
            self.fc3 = nn.Linear(8, 1)
            self.sigmoid = nn.Sigmoid()

        def forward(self, x):
            x = self.relu1(self.fc1(x))
            x = self.relu2(self.fc2(x))
            x = self.sigmoid(self.fc3(x))
            return x
else:
    DeepHistovarNet = None

class MLScorer:
    def __init__(self, model_type="heuristic", model_path=None):
        self.model_type = model_type
        self.is_trained = False
        self.feature_names = ["depth", "af", "mapq", "strand_bias", "pos_bias", "posterior_gap"]
        self.model = None
        
        if model_type == "heuristic":
            pass
        elif model_type == "rf":
            if not HAS_SKLEARN:
                raise ImportError("scikit-learn is required for model_type='rf'")
            self.model = RandomForestClassifier(n_estimators=100, random_state=42)
        elif model_type == "boosting":
            if not HAS_SKLEARN:
                raise ImportError("scikit-learn is required for model_type='boosting'")
            self.model = GradientBoostingClassifier(n_estimators=100, random_state=42)
        elif model_type == "lr":
            if not HAS_SKLEARN:
                raise ImportError("scikit-learn is required for model_type='lr'")
            self.model = LogisticRegression(random_state=42)
        elif model_type == "dl":
            if not HAS_TORCH:
                raise ImportError("PyTorch is required for model_type='dl'")
            self.model = DeepHistovarNet(len(self.feature_names))
            self.optimizer = optim.Adam(self.model.parameters(), lr=0.01)
            self.criterion = nn.BCELoss()
        else:
            raise ValueError(f"Unknown model type: {model_type}")

        if model_path:
            self.load(model_path)

    def save(self, path):
        if self.model_type == "heuristic":
            raise ValueError("Heuristic scorer has no model weights to save")
        if self.model_type == "dl":
            torch.save(self.model.state_dict(), path)
        else:
            with open(path, "wb") as f:
                pickle.dump(self.model, f)

    def load(self, path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Model file not found: {path}")
        if self.model_type == "dl":
            state = torch.load(path, map_location="cpu", weights_only=True)
            if isinstance(state, dict):
                self.model.load_state_dict(state)
            else:
                self.model = state
            self.model.eval()
        elif self.model_type in ("rf", "boosting", "lr"):
            with open(path, "rb") as f:
                self.model = pickle.load(f)
        self.is_trained = True

    def _dict_to_array(self, features_list):
        X = []
        for f in features_list:
            x_row = [f.get(k, 0) for k in self.feature_names]
            X.append(x_row)
        return np.array(X)

    def train(self, X_dicts, y):
        if self.model_type == "heuristic":
            raise ValueError("Cannot train model_type='heuristic'")
        X = self._dict_to_array(X_dicts)
        
        if self.model_type == "dl":
            # PyTorch training loop
            X_tensor = torch.tensor(X, dtype=torch.float32)
            y_tensor = torch.tensor(y, dtype=torch.float32).unsqueeze(1)
            
            self.model.train()
            for epoch in range(50):
                self.optimizer.zero_grad()
                outputs = self.model(X_tensor)
                loss = self.criterion(outputs, y_tensor)
                loss.backward()
                self.optimizer.step()
        else:
            self.model.fit(X, y)
            
        self.is_trained = True
        
    def predict_proba(self, X_dicts):
        if not self.is_trained:
            # Fallback heuristic probability for uninitialized/mock runs instead of returning exactly 0.98
            probas = []
            for d in X_dicts:
                af = d.get("af", 0)
                depth = d.get("depth", 0)
                mapq = d.get("mapq", 0)
                post_gap = d.get("posterior_gap", 1.0)
                
                # Continuous heuristic score
                # If uniquely mapped (post_gap >= 0.2), we expect high MAPQ.
                # If multi-mapped (post_gap < 0.2), MAPQ is expected to be low/0, so we do not penalize it.
                if post_gap < 0.2:
                    score = af * min(depth / 20.0, 1.0)
                else:
                    score = af * min(depth / 20.0, 1.0) * min(mapq / 40.0, 1.0)
                # Scale it so a good variant gets > 0.5
                score = min(0.99, max(0.01, score * 5.0))
                
                probas.append(score)
            return probas
            
        X = self._dict_to_array(X_dicts)
        
        if self.model_type == "dl":
            self.model.eval()
            with torch.no_grad():
                X_tensor = torch.tensor(X, dtype=torch.float32)
                outputs = self.model(X_tensor)
                return outputs.squeeze(1).tolist()
        else:
            probas = self.model.predict_proba(X)
            # Return probability of class 1 (true variant)
            return probas[:, 1].tolist()

def score():
    # Keep original dummy function
    return 0.98
