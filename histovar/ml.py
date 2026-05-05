from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim

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

class MLScorer:
    def __init__(self, model_type="rf"):
        self.model_type = model_type
        self.is_trained = False
        self.feature_names = ["depth", "af", "mapq", "strand_bias", "pos_bias", "posterior_gap"]
        
        if model_type == "rf":
            self.model = RandomForestClassifier(n_estimators=100, random_state=42)
        elif model_type == "boosting":
            self.model = GradientBoostingClassifier(n_estimators=100, random_state=42)
        elif model_type == "lr":
            self.model = LogisticRegression(random_state=42)
        elif model_type == "dl":
            self.model = DeepHistovarNet(len(self.feature_names))
            self.optimizer = optim.Adam(self.model.parameters(), lr=0.01)
            self.criterion = nn.BCELoss()
        else:
            raise ValueError("Unknown model type")

    def _dict_to_array(self, features_list):
        X = []
        for f in features_list:
            x_row = [f.get(k, 0) for k in self.feature_names]
            X.append(x_row)
        return np.array(X)

    def train(self, X_dicts, y):
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
                if depth >= 10 and af >= 0.05:
                    probas.append(0.95)
                elif depth >= 5 and af >= 0.01:
                    probas.append(0.80)
                else:
                    probas.append(0.01)
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
