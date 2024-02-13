import numpy as np
import pandas as pd


def CR2_variance_estimator(model, clustid):
    """
    Calculates the CR2 cluster-robust variance estimator.

    Args:
        model: A fitted statsmodels model object.
        clustid: The name of the column in the model's DataFrame that 
                 identifies the clusters.

    Returns:
        A pandas Series containing the CR2 robust standard errors.
    """

    df = model.model.data.orig_endog.index.to_frame().join(model.model.data.frame)
    cluster_groups = df.groupby(clustid)

    # Matrices from the equations 
    X = model.model.exog  
    y = model.model.endog
    W = np.diag(cluster_groups.size().values) 
    Phi = _get_cluster_error_covariance(model, cluster_groups)

    # Terms from the CR2 adjustments
    Z = model.model.wexog 
    Mz = np.linalg.pinv(Z.T @ W @ Z)
    Hz = Z @ Mz @ Z.T @ W
    A = []
    B = []

    for i, (cluster_id, group) in enumerate(cluster_groups):
        Ci = _get_selection_matrix(group.index.values, X.shape[0])
        Di = np.linalg.cholesky(Phi[i]) 
        S_bar = Ci @ X - Ci @ Hz @ X
        U_tilde = S_bar @ np.delete(model.params, np.argwhere(np.isnan(model.params)), axis=0)
        R_dot = U_tilde[:, 0, None] 

        A.append(R_dot.T @ np.linalg.inv(Di.T @ Di) @ R_dot)
        B.append(Di.T @ Ci @ W @ Ci.T @ Di)

    A = np.concatenate(A).squeeze()
    B = np.linalg.inv(np.block_diag(*B)) 

    #  CR2 estimator calculation
    e = model.resid
    V_CR2 = Mz @ Z.T @ W @ np.diag(e**2) @ W @ Z @ Mz 
    V_CR2 += Mz @ Z.T @ W @ A @ W @ Z @ Mz
    V_CR2 += Mz @ Z.T @ W @ B @ W @ Z @ Mz

    robust_se = np.sqrt(np.diag(V_CR2)) 
    return pd.Series(robust_se, index=model.params.index)

# Helper functions
def _get_selection_matrix(indices, n_obs):
    C = np.zeros((len(indices), n_obs))
    C[np.arange(len(indices)), indices] = 1
    return C

def _get_cluster_error_covariance(model, cluster_groups):
    Phi = []
    for _, group in cluster_groups:
        e = model.resid[group.index]
        Phi.append(e[:, None] @ e[None, :])
    return np.asarray(Phi) 
