#!/usr/bin/env python3
"""
Master Module 4: Gated Attention Deep Multiple Instance Learning (Attention-MIL).
Translates Whole-Slide Diagnostic H&E morphology (< 5 €) into microvascular decoupling
and progression-free survival risk scores without spatial transcriptomics.
"""

import math
from typing import Dict, List, Optional, Tuple
import numpy as np


class GatedAttentionLayer:
    """Gated Attention Pooling layer for bag-of-features whole-slide representations.
    
    Formula:
        a_k = softmax(w^T (tanh(V h_k) * sigmoid(U h_k)))
    """
    def __init__(self, in_features: int = 512, hidden_dim: int = 128, seed: int = 42):
        rng = np.random.RandomState(seed)
        scale_v = math.sqrt(2.0 / (in_features + hidden_dim))
        scale_u = math.sqrt(2.0 / (in_features + hidden_dim))
        scale_w = math.sqrt(2.0 / (hidden_dim + 1))
        
        self.V = rng.normal(0, scale_v, size=(hidden_dim, in_features))
        self.b_v = np.zeros(hidden_dim)
        self.U = rng.normal(0, scale_u, size=(hidden_dim, in_features))
        self.b_u = np.zeros(hidden_dim)
        self.w = rng.normal(0, scale_w, size=(1, hidden_dim))
        self.b_w = np.zeros(1)
        
    def forward(self, H: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Computes patch attention weights and slide-level aggregated embedding.
        
        Args:
            H: Patch feature matrix of shape (N_patches, in_features).
            
        Returns:
            z: Slide representation vector of shape (in_features,).
            attention_weights: Array of shape (N_patches,) summing to 1.0.
        """
        # Linear projections
        v_out = np.tanh(np.dot(H, self.V.T) + self.b_v)
        u_out = 1.0 / (1.0 + np.exp(-np.clip(np.dot(H, self.U.T) + self.b_u, -15, 15)))
        
        # Gating mechanism
        gated = v_out * u_out
        scores = np.dot(gated, self.w.T).squeeze(-1) + self.b_w
        
        # Softmax normalization
        scores_exp = np.exp(scores - np.max(scores))
        attention_weights = scores_exp / np.maximum(np.sum(scores_exp), 1e-12)
        
        z = np.dot(attention_weights, H)
        return z, attention_weights


class AttentionSurvivalMIL:
    """Deep Multiple Instance Learning architecture for survival prediction."""
    def __init__(self, in_features: int = 512, hidden_dim: int = 128, seed: int = 42):
        self.attention = GatedAttentionLayer(in_features, hidden_dim, seed=seed)
        rng = np.random.RandomState(seed + 1)
        self.classifier_w = rng.normal(0, 0.05, size=in_features)
        self.classifier_b = 0.0
        
    def predict_hazard(self, bag_features: np.ndarray) -> Tuple[float, np.ndarray]:
        """Predicts log-hazard ratio theta and patch attention weights."""
        z, attention = self.attention.forward(bag_features)
        log_hazard = float(np.dot(z, self.classifier_w) + self.classifier_b)
        return log_hazard, attention
    
    def predict_cohort(self, bags: List[np.ndarray]) -> np.ndarray:
        """Predicts log-hazard ratios across a cohort of patient bags."""
        hazards = np.zeros(len(bags))
        for i, bag in enumerate(bags):
            hazards[i], _ = self.predict_hazard(bag)
        return hazards


def cox_negative_log_likelihood(hazards: np.ndarray, times: np.ndarray, events: np.ndarray) -> float:
    """Computes Cox Proportional Hazards Negative Partial Log-Likelihood."""
    order = np.argsort(-times)
    sorted_hazards = hazards[order]
    sorted_events = events[order]
    
    risk_sum = 0.0
    loss = 0.0
    for i in range(len(sorted_hazards)):
        risk_sum += np.exp(sorted_hazards[i])
        if sorted_events[i] == 1:
            loss += np.log(max(risk_sum, 1e-12)) - sorted_hazards[i]
            
    return float(loss / max(np.sum(events), 1))


def project_spatial_attention_heatmap(
    coords: np.ndarray,
    attention_weights: np.ndarray,
    grid_size: int = 100,
    sigma: float = 1.5
) -> np.ndarray:
    """Projects 1D patch attention weights onto a 2D spatial tissue grid with Gaussian smoothing."""
    from scipy.ndimage import gaussian_filter
    
    x_min, x_max = np.min(coords[:, 0]), np.max(coords[:, 0])
    y_min, y_max = np.min(coords[:, 1]), np.max(coords[:, 1])
    
    x_norm = np.clip(((coords[:, 0] - x_min) / max(x_max - x_min, 1e-6) * (grid_size - 1)).astype(int), 0, grid_size - 1)
    y_norm = np.clip(((coords[:, 1] - y_min) / max(y_max - y_min, 1e-6) * (grid_size - 1)).astype(int), 0, grid_size - 1)
    
    grid = np.zeros((grid_size, grid_size), dtype=float)
    scaled_attn = attention_weights / max(np.max(attention_weights), 1e-12)
    for i in range(len(coords)):
        grid[y_norm[i], x_norm[i]] += scaled_attn[i]
        
    return gaussian_filter(grid, sigma=sigma)
