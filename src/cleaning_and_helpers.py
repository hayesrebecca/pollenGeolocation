import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error

def clean_rbcl_data(data_filepath,
                    save_filepath):
    df = pd.read_csv(data_filepath, low_memory=False)
    df['GenSp'] = df['Genus'] + '.' + df['Species']
    rbcl_cols = [col for col in df if col.startswith('RBCL')]
    my_cols = ['UniqueID', 'Family', 'Genus', 'GenSp', 'Site', 'Lat', 'Long'] + rbcl_cols
    df = df.reindex(columns=my_cols)
    df = df.dropna(subset=[i for i in rbcl_cols], how='all')
    df = df.fillna(0)
    df.to_csv(save_filepath)
    return df

def plot_test_preds(y_test, preds, scaler, model_type, ax):
    # Back-transform
    unscaled_y_test = scaler.inverse_transform(y_test)
    unscaled_preds = scaler.inverse_transform(preds)

    y_test_split = np.hsplit(unscaled_y_test, 2)
    preds_split = np.hsplit(unscaled_preds, 2)

    # Calculate R^2 and MSE
    r2 = r2_score(y_test, preds)
    mse = mean_squared_error(y_test, preds)

    # Scatter plots
    ax.scatter(y_test_split[0], y_test_split[1], marker='*', s=200, label='Real Location', color='orange')
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
