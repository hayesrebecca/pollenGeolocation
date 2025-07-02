import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error

def clean_rbcl_data(data_filepath,
                    save_filepath,
                    format="utf-8"):
    df = pd.read_csv(data_filepath, low_memory=False, encoding=format)

    if "UniqueID" in df.columns:
        df['GenSp'] = df['Genus'] + '.' + df['Species']
        rbcl_cols = [col for col in df if col.startswith('RBCL')]
        my_cols = ['UniqueID', 'Family', 'Genus', 'GenSp', 'Site', 'Lat', 'Long'] + rbcl_cols
        df = df.reindex(columns=my_cols)
        df = df.dropna(subset=[i for i in rbcl_cols], how='all')
        df = df.fillna(0)
        df.to_csv(save_filepath)
    else:
        df["UniqueID"] = df["SampleID"]
        df['GenSp'] = df['Genus'] + '.' + df['Species']
        rbcl_cols = [col for col in df if col.startswith('RBCL')]
        my_cols = ['UniqueID', 'Family', 'Genus', 'GenSp', 'Site', 'Lat', 'Long'] + rbcl_cols
        df = df.reindex(columns=my_cols)
        df = df.dropna(subset=[i for i in rbcl_cols], how='all')
        df = df.fillna(0)
        df.to_csv(save_filepath)
    return df

""" def plot_test_preds(y_test, preds, scaler, model_type, ax):
    # Back-transform
    unscaled_y_test = scaler.inverse_transform(y_test)
    unscaled_preds = scaler.inverse_transform(preds)

    y_test_split = np.hsplit(unscaled_y_test, 2)
    preds_split = np.hsplit(unscaled_preds, 2)

    # Calculate R^2 and MSE
    r2 = r2_score(y_test, preds)
    mse = mean_squared_error(y_test, preds)

    # Scatter plots
    ax.scatter(y_test_split[0], y_test_split[1], marker='*', s=200, label='Real Location', color='green')
    ax.scatter(preds_split[0], preds_split[1], alpha=0.8, label='Predicted Location', color='blue')

    # Customize to match `theme_linedraw()`
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.tick_params(axis='both', which='major', labelsize=10, color='black')
    ax.grid(True, color='lightgray', linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)

    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)

    # Add labels
    ax.set_xlabel('Latitude', fontsize=12, color='black') 
    ax.set_ylabel('Longitude', fontsize=12, color='black') 

    # Add legend
    ax.legend(frameon=True, fontsize=10)

    # Add text for R^2 and MSE
    ax.text(
        0.05, 0.20,
        f'{model_type}\n$R^2$: {r2:.2f}\nMSE: {mse:.2f}',
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.3', edgecolor='gray', facecolor='white')
    )

    ax.invert_xaxis()
    ax.invert_yaxis() """

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, mean_squared_error
from geopy.distance import geodesic
import matplotlib.colors as mcolors
from string import ascii_uppercase

def plot_test_preds(y_test, preds, scaler, model_type, ax, norm):
    # Back-transform
    unscaled_y_test = scaler.inverse_transform(y_test)
    unscaled_preds = scaler.inverse_transform(preds)

    # Extract latitude and longitude separately
    lat_test, lon_test = unscaled_y_test[:, 0], unscaled_y_test[:, 1]
    lat_pred, lon_pred = unscaled_preds[:, 0], unscaled_preds[:, 1]

    # Calculate great-circle distance (in km) for each point
    distances = np.array([
        geodesic((lat1, lon1), (lat2, lon2)).kilometers
        for lat1, lon1, lat2, lon2 in zip(lat_test, lon_test, lat_pred, lon_pred)
    ])

    # Sort points by distance (smallest first) for proper layering
    sorted_indices = np.argsort(distances)
    lat_pred, lon_pred = lat_pred[sorted_indices], lon_pred[sorted_indices]
    distances = distances[sorted_indices]

    # Calculate R^2 and MSE
    r2 = r2_score(y_test, preds)
    mse = mean_squared_error(y_test, preds)

    # Scatter plots with black edge for visibility
    ax.scatter(lat_test, lon_test, marker='*', s=350, label='Real Location', 
               color='gold', edgecolors='black', linewidths=0.8)
    sc = ax.scatter(lat_pred, lon_pred, c=distances, cmap='magma', alpha=0.9, 
                    label='Predicted Location', edgecolors='black', linewidths=0.6, norm=norm)

    # Customize to match `theme_linedraw()`
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.tick_params(axis='both', which='major', labelsize=10, color='black')
    ax.grid(True, color='lightgray', linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)

    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)

    # Add labels
    ax.set_xlabel('Latitude', fontsize=12, color='black') 
    ax.set_ylabel('Longitude', fontsize=12, color='black') 

    # Add legend
    ax.legend(frameon=True, fontsize=10)

    # Add text for R^2 and MSE
    ax.text(
        0.02, 0.25,
        f'{model_type}\n$R^2$: {r2:.2f}\nMSE: {mse:.2f}',
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.3', edgecolor='gray', facecolor='white')
    )

    ax.invert_xaxis()
    ax.invert_yaxis()

    return sc  # Return scatter plot for colorbar reference

# -----------------------------
# Function to split each project
# -----------------------------
def split_project(X, y, test_size, random_state):
    return train_test_split(X, y, test_size=test_size, random_state=random_state)

# -----------------------------
# Custom RMSE scorer
# -----------------------------
def multioutput_rmse(y_true, y_pred):
    return np.sqrt(((y_true - y_pred) ** 2).mean())

from sklearn.metrics import r2_score, mean_squared_error, median_absolute_error
from geopy.distance import geodesic
import numpy as np

def evaluate_model(name, model_class, best_params, X_train, y_train, X_test, y_test):
    """
    Initialize, fit, predict, and evaluate a model using best hyperparameters.

    Parameters:
    - name: string, name of model
    - model_class: class or callable that returns a model
    - best_params: dict of hyperparameters
    - X_train, y_train, X_test, y_test: data splits

    Returns:
    - dict with performance metrics and predictions
    """
    print(f"Evaluating {name}...")

    # If model_class is a callable (e.g., lambda), call it to get the estimator
    model = model_class(**best_params)

    model.fit(X_train, y_train)
    preds = model.predict(X_test)

    r2 = r2_score(y_test, preds)
    rmse = mean_squared_error(y_test, preds, squared=False)
    mae = median_absolute_error(y_test, preds)
    distances = [geodesic(real, pred).kilometers for real, pred in zip(y_test, preds)]
    avg_dist = np.mean(distances)
    se_dist = np.std(distances)

    return {
        "model": name,
        "r2": r2,
        "rmse": rmse,
        "mae": mae,
        "avg_km_error": avg_dist,
        "se_km_error": se_dist,
        "preds": preds
    }

from sklearn.metrics import r2_score, mean_squared_error, median_absolute_error
from geopy.distance import geodesic
import numpy as np
import pandas as pd

def evaluate_model_per_project(name, model_class, X_train, y_train, X_test, y_test, project_test=None):
    """
    Initialize, fit, predict, and evaluate a model with default hyperparameters per project.
    """
    print(f"Evaluating {name}...")

    # If it’s a lambda or a class, calling it gets you a fresh model
    model = model_class()

    model.fit(X_train, y_train)
    preds = model.predict(X_test)

    # Overall metrics
    r2 = r2_score(y_test, preds)
    rmse = mean_squared_error(y_test, preds, squared=False)
    mae = median_absolute_error(y_test, preds)
    distances = [geodesic(real, pred).kilometers for real, pred in zip(y_test, preds)]
    avg_dist = np.mean(distances)
    se_dist = np.std(distances)

    results = {
        "model": name,
        "r2": r2,
        "rmse": rmse,
        "mae": mae,
        "avg_km_error": avg_dist,
        "se_km_error": se_dist,
        "preds": preds
    }

    # ---- Compute per-project metrics ----
    if project_test is not None:
        df = pd.DataFrame({
            "y_true_lat": [y[0] for y in y_test],
            "y_true_lon": [y[1] for y in y_test],
            "y_pred_lat": [y[0] for y in preds],
            "y_pred_lon": [y[1] for y in preds],
            "project": project_test
        })

        project_metrics = []

        for project_id, group in df.groupby("project"):
            y_true_proj = group[["y_true_lat", "y_true_lon"]].values
            y_pred_proj = group[["y_pred_lat", "y_pred_lon"]].values

            r2_p = r2_score(y_true_proj, y_pred_proj)
            rmse_p = mean_squared_error(y_true_proj, y_pred_proj, squared=False)
            mae_p = median_absolute_error(y_true_proj, y_pred_proj)
            distances_p = [geodesic(real, pred).kilometers for real, pred in zip(y_true_proj, y_pred_proj)]
            avg_dist_p = np.mean(distances_p)
            se_dist_p = np.std(distances_p)

            project_metrics.append({
                "project": project_id,
                "r2": r2_p,
                "rmse": rmse_p,
                "mae": mae_p,
                "avg_km_error": avg_dist_p,
                "se_km_error": se_dist_p
            })

        results["per_project"] = project_metrics

    return results


