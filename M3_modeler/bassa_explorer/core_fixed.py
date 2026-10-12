"""Test copy of bassa_reg/spike_and_slab/spike_and_slab_core.py (github.com/Milo-group/BASSA, commit 661a164): PATCH 1, 2 and 3.
The package itself is not modified. The released core is wrong when descriptors are correlated: reconsidering a descriptor
that is in the model, it reads the determinant ratio off the wrong row of the precision matrix; this copy matches the exact
posterior (sampler_check.py in the DescriPyTor paper's notes/analysis/bassa)."""
import numpy as np

from numba import float32, int32, njit
from numba.experimental import jitclass

# import numba
# numba.config.DISABLE_JIT = True
# disable numba

@njit(fastmath=True)
def build_xtx_xty_pairs(x, y, s0, new_feat_idx):
    x = np.ascontiguousarray(x)
    y = np.ascontiguousarray(y)

    n_rows, n_feats = x.shape
    k = 0
    for j in range(n_feats):
        if s0[j]:
            k += 1

    xs0 = np.zeros((n_rows, k), dtype=np.float32)
    col_idx = 0
    for j in range(n_feats):
        if s0[j]:
            xs0[:, col_idx] = x[:, j]
            col_idx += 1
    xs0 = np.ascontiguousarray(xs0)

    xtx_0 = xs0.T @ xs0
    xty_0 = xs0.T @ y

    v = np.ascontiguousarray(x[:, int(new_feat_idx)])
    vtv = v @ v
    vtx = v.T @ xs0
    vty = v @ y

    xtx_1 = np.zeros((k + 1, k + 1), dtype=np.float32)
    xtx_1[0, 0] = vtv
    xtx_1[0, 1:] = vtx
    xtx_1[1:, 0] = vtx
    xtx_1[1:, 1:] = xtx_0

    xty_1 = np.zeros(k + 1, dtype=np.float32)
    xty_1[0] = vty
    xty_1[1:] = xty_0

    return xtx_0, xtx_1, xty_0, xty_1

@njit(fastmath=True)
def build_xtx_xty_pair(x, y, s0):
    x = np.ascontiguousarray(x)
    y = np.ascontiguousarray(y)

    n_rows, n_feats = x.shape
    k = 0
    for j in range(n_feats):
        if s0[j]:
            k += 1

    xs0 = np.zeros((n_rows, k), dtype=np.float32)
    col_idx = 0
    for j in range(n_feats):
        if s0[j]:
            xs0[:, col_idx] = x[:, j]
            col_idx += 1
    xs0 = np.ascontiguousarray(xs0)

    xtx_0 = xs0.T @ xs0
    xty_0 = xs0.T @ y

    return xtx_0, xty_0

@njit(fastmath=True)
def move_front_back(xtx, xty, dest_pos):
    k = xtx.shape[0]
    if dest_pos == 0:                       #
        return xtx, xty

    perm = np.empty(k, dtype=np.int32)
    # first, rows/cols 1 … k-1
    for i in range(dest_pos):
        perm[i] = i + 1
    perm[dest_pos] = 0
    for i in range(dest_pos + 1, k):
        perm[i] = i

    xtx_new = np.empty_like(xtx)
    xty_new = np.empty_like(xty)
    for i in range(k):
        src_i = perm[i]
        xty_new[i] = xty[src_i]
        for j in range(k):
            xtx_new[i, j] = xtx[src_i, perm[j]]

    return xtx_new, xty_new

@njit(fastmath=True)
def build_xtx_xty_pair_for_extra_feature(x, y, xtx_0, xty_0, s0, new_feat_idx):
    x = np.ascontiguousarray(x)
    y = np.ascontiguousarray(y)

    n_rows, n_feats = x.shape
    k = 0
    for j in range(n_feats):
        if s0[j]:
            k += 1

    xs0 = np.zeros((n_rows, k), dtype=np.float32)
    col_idx = 0
    for j in range(n_feats):
        if s0[j]:
            xs0[:, col_idx] = x[:, j]
            col_idx += 1
    xs0 = np.ascontiguousarray(xs0)

    v = np.ascontiguousarray(x[:, int(new_feat_idx)])
    vtv = v @ v
    vtx = v.T @ xs0
    vty = v @ y

    xtx_1 = np.zeros((k + 1, k + 1), dtype=np.float32)
    xtx_1[0, 0] = vtv
    xtx_1[0, 1:] = vtx
    xtx_1[1:, 0] = vtx
    xtx_1[1:, 1:] = xtx_0

    xty_1 = np.zeros(k + 1, dtype=np.float32)
    xty_1[0] = vty
    xty_1[1:] = xty_0

    return xtx_1, xty_1

@njit(fastmath=True)
def get_values_for_s(s, xs, y, tau, sigma, xtx, xty):
    s_sum = int(np.sum(s))
    if s_sum == 0:
        return np.zeros(1, dtype=np.float32), np.empty((0, 0), dtype=np.float32)

    if tau == 0:
        sigma_inverse = (sigma ** -2) * xtx
    else:
        sigma_inverse = (tau ** -2) * np.eye(int(np.sum(s)), dtype=np.float32) + (sigma ** -2) * xtx

    V = np.ascontiguousarray(xty.T)
    chol_sol = solve_with_cholesky(sigma_inverse, V)
    ytx_sig_xty = np.ascontiguousarray(chol_sol @ xty)

    return ytx_sig_xty, sigma_inverse

@njit(fastmath=True)
def solve_with_cholesky(A, V):
    epsilon = np.float32(1e-6)
    A = A + epsilon * np.eye(A.shape[0], dtype=np.float32)
    L = np.linalg.cholesky(A)
    y = np.linalg.solve(L, V)
    x = np.linalg.solve(L.T, y)
    return x

@njit(fastmath=True)
def get_determinants_ratio(sigma1, sigma0):
    a11 = sigma1[0][0]
    v = sigma1[0][1:]
    L = np.linalg.cholesky(sigma0)
    y = np.linalg.solve(L, v.T)
    vtsig = y.T @ y
    new_ratio = 1 / (a11 - vtsig)

    return new_ratio

@njit(fastmath=True)
def create_feature_matrices(x, s0, new_feature_index):
    n_rows = x.shape[0]
    n_selected = 0
    for i in range(s0.shape[0]):
        if s0[i]:
            n_selected += 1

    xs_0 = np.zeros((n_rows, n_selected), dtype=np.float32)
    idx = 0
    for i in range(s0.shape[0]):
        if s0[i]:
            xs_0[:, idx] = x[:, i]
            idx += 1

    xs_1 = np.zeros((n_rows, n_selected + 1), dtype=np.float32)
    xs_1[:, 0] = x[:, int(new_feature_index)]
    for j in range(n_selected):
        xs_1[:, j + 1] = xs_0[:, j]

    return np.ascontiguousarray(xs_0), np.ascontiguousarray(xs_1)

@njit(fastmath=True)
def get_xs(s, x):
    indices = np.nonzero(s)[0]
    xs = x[:, indices]
    return xs

@njit(fastmath=True)
def insert_masked(mask, w):
    output = np.zeros(len(mask), dtype=np.float32)
    w_index = 0
    for i in range(len(mask)):
        if mask[i] == 1:
            output[i] = w[w_index]
            w_index += 1
    return output

spec = [
    # bassa_reg state
    ('x', float32[:, :]),
    ('y', float32[:]),
    ('x_test', float32[:, :]),
    ('y_test', float32[:]),
    ('s', float32[:]),
    ('w', float32[:]),
    ('a', float32),
    ('sigma', float32),
    ('tau', float32),
    ('a_alpha_prior', float32),
    ('a_beta_prior', float32),
    ('sigma_alpha_prior', float32),
    ('sigma_beta_prior', float32),
    ('tau_alpha_prior', float32),
    ('tau_beta_prior', float32),

    # bookkeeping
    ('t', int32),
    ('num_iters', int32),

    # history (with dynamic shape)
    ('s_chain', float32[:, :]),
    ('a_chain', float32[:]),
    ('sigma_chain', float32[:]),
    ('tau_chain', float32[:]),
    ('w_chain', float32[:, :]),
    ('w_dim_chain', int32[:]),
    ('internal_y_predictions', float32[:, :]),
    ('y_predictions', float32[:, :]),

    # helpers
    ('exp_subtraction_threshold', float32),

    # Geweke's test
    ('geweke_sigma', float32[:, :]),
    ('geweke_tau', float32[:, :]),
    ('geweke_a', float32[:, :]),
    ('geweke_s_sum', float32[:, :]),
    ('geweke_index', int32),
]

@jitclass(spec)
class SpikeAndSlabCore:
    def __init__(self, x, y, x_test, class_priors,
                 num_iters, prediction_iters):
        self.x = np.ascontiguousarray(x).astype(np.float32)
        self.y = np.ascontiguousarray(y).astype(np.float32)
        self.x_test = np.ascontiguousarray(x_test).astype(np.float32)

        max_w_dim = x.shape[1]
        self.num_iters = num_iters
        self.t = 0

        self.a_alpha_prior, self.a_beta_prior, \
        self.sigma_alpha_prior, self.sigma_beta_prior, \
        self.tau_alpha_prior, self.tau_beta_prior = class_priors

        n_features = x.shape[1]
        self.s = np.zeros(n_features, dtype=np.float32)
        self.s = np.ascontiguousarray(self.s).astype(np.float32)
        self.w = np.zeros(max_w_dim, dtype=np.float32)
        self.w = np.ascontiguousarray(self.w).astype(np.float32)
        self.a = np.float32(0.0)
        self.sigma = np.float32(1.0)
        self.tau = np.float32(1.0)
        self.exp_subtraction_threshold = np.float32(50)

        # Initialize predictions
        n_points = self.x.shape[0]
        self.internal_y_predictions = np.zeros((n_points, num_iters + 1), dtype=np.float32)
        self.y_predictions = np.zeros((self.x_test.shape[0], prediction_iters), dtype=np.float32)

        #   Initialize chain arrays
        self.s_chain = np.zeros((n_features, num_iters), dtype=np.float32)
        self.a_chain = np.zeros(num_iters, dtype=np.float32)
        self.sigma_chain = np.zeros(num_iters, dtype=np.float32)
        self.tau_chain = np.zeros(num_iters, dtype=np.float32)
        self.w_chain = np.zeros((max_w_dim, num_iters), dtype=np.float32)
        self.w_dim_chain = np.zeros(num_iters, dtype=np.int32)


        #   Initialize Geweke's test arrays
        self.geweke_sigma = np.zeros((num_iters, 2), dtype=np.float32)
        self.geweke_tau = np.zeros((num_iters, 2), dtype=np.float32)
        self.geweke_a = np.zeros((num_iters, 2), dtype=np.float32)
        self.geweke_s_sum = np.zeros((num_iters, 2), dtype=np.float32)
        self.geweke_index = 0

    def set_initial_state(self, a_val, s_val, w_val, tau_val, sigma_val):
        self.a = np.float32(a_val)
        self.s = np.ascontiguousarray(s_val).astype(np.float32)
        self.tau = np.float32(tau_val)
        self.sigma = np.float32(sigma_val)
        self.w = np.ascontiguousarray(w_val).astype(np.float32)

    def seed_random(self, seed):
        np.random.seed(seed)

    def predict_on_test_set(self, iteration_idx, samples_per_y):
        xs_test = get_xs(self.s, self.x_test)
        s_size = int(np.sum(self.s))

        if s_size == 0:
            y_pred_mean = np.zeros(self.x_test.shape[0], dtype=np.float32)
        else:
            w_used = self.w[:xs_test.shape[1]]
            w_used = np.ascontiguousarray(w_used)
            y_pred_mean = (xs_test @ w_used).astype(np.float32)

        self.y_predictions[:, iteration_idx] = y_pred_mean

    def compute_internal_y_prediction(self):
        xs = get_xs(self.s, self.x)
        s_size = int(np.sum(self.s))

        if s_size == 0:
            y_pred_mean = np.zeros(self.x.shape[0], dtype=np.float32)
        else:
            w_used = self.w
            w_used = np.ascontiguousarray(w_used)
            y_pred_mean = (xs @ w_used).astype(np.float32)

        self.internal_y_predictions[:, self.t] = y_pred_mean

    def sample_geweke(self):
        sigma0 = np.float32(np.random.gamma(self.sigma_alpha_prior, scale=1.0 / self.sigma_beta_prior) ** -0.5)
        tau0 = np.float32(np.random.gamma(self.tau_alpha_prior, scale=1.0 / self.tau_beta_prior) ** -0.5)
        a0 = np.float32(np.random.beta(self.a_alpha_prior, self.a_beta_prior))
        s0 = np.random.binomial(1, a0, size=self.x.shape[1]).astype(np.float32)
        w0 = (np.random.normal(0.0, tau0, size=self.x.shape[1]) * s0).astype(np.float32)

        x_c = np.ascontiguousarray(self.x)
        noise = np.random.normal(0.0, sigma0, size=self.x.shape[0]).astype(np.float32)
        y0 = (x_c @ w0 + noise).astype(np.float32)

        _y = self.y.copy()
        _sigma = self.sigma
        _tau = self.tau
        _a = self.a
        _s = self.s.copy()
        _w = self.w.copy()

        self.y = y0
        self.sigma = sigma0
        self.tau = tau0
        self.a = a0
        self.s = s0
        w_compact = np.zeros(int(np.sum(s0)), dtype=np.float32)
        k = 0
        for i in range(s0.shape[0]):
            if s0[i] == 1.0:
                w_compact[k] = w0[i]
                k += 1
        self.w = w_compact

        self.cycle(False)

        sigma1 = self.sigma
        tau1 = self.tau
        a1 = self.a
        s1_sum = np.sum(self.s)

        self.y = _y
        self.sigma = _sigma
        self.tau = _tau
        self.a = _a
        self.s = _s
        self.w = _w

        i = self.geweke_index

        self.geweke_sigma[i, 0] = sigma0
        self.geweke_sigma[i, 1] = sigma1

        self.geweke_tau[i, 0] = tau0
        self.geweke_tau[i, 1] = tau1

        self.geweke_a[i, 0] = a0
        self.geweke_a[i, 1] = a1

        self.geweke_s_sum[i, 0] = np.sum(s0)
        self.geweke_s_sum[i, 1] = s1_sum

        self.geweke_index += 1

    def sample_s(self, s_index_to_sample, xtx_0, xty_0, xtx_1, xty_1, inserted_0):
        s0 = self.s.copy()
        s0[int(s_index_to_sample)] = np.float32(0)
        s1 = self.s.copy()
        s1[int(s_index_to_sample)] = np.float32(1)

        if inserted_0:
            xtx_1, xty_1 = build_xtx_xty_pair_for_extra_feature(self.x, self.y, xtx_0, xty_0, s0, s_index_to_sample)
            xtx_1f, xty_1f = xtx_1, xty_1
        else:
            xtx_0, xty_0 = build_xtx_xty_pair(self.x, self.y, s0)
            # PATCH 3: get_determinants_ratio() reads row 0 of the larger precision matrix as the candidate descriptor.
            # When the descriptor is already in the model the matrix handed in is in natural order, so row 0 is the
            # first active descriptor, not the candidate. Rebuild the larger matrix with the candidate in front.
            xtx_1f, xty_1f = build_xtx_xty_pair_for_extra_feature(self.x, self.y, xtx_0, xty_0, s0, s_index_to_sample)

        xs_0, xs_1 = create_feature_matrices(self.x, s0, s_index_to_sample)
        ytx_sig_xty_0, sigma_inverse_0 = get_values_for_s(s0, xs_0, self.y, self.tau, self.sigma, xtx_0, xty_0)
        ytx_sig_xty_1, sigma_inverse_1 = get_values_for_s(s1, xs_1, self.y, self.tau, self.sigma, xtx_1f, xty_1f)
        det_ratio = np.float32(0.0)
        exp = np.float32(0.0)
        pre_exp_expression = np.float32(0.0)
        exp_subtraction = (ytx_sig_xty_1 - ytx_sig_xty_0) / (np.float32(2) * (self.sigma ** np.float32(4)))
        s1_prob = np.float32(0.0)
        if self.exp_subtraction_threshold >= exp_subtraction >= -self.exp_subtraction_threshold:
            det_ratio = get_determinants_ratio(sigma_inverse_1, sigma_inverse_0)
            det_ratio = det_ratio ** np.float32(0.5)
            det_ratio = det_ratio / self.tau
            # pre_exp_expression = det_ratio * self.a / ((1 - self.a) * self.tau)
            pre_exp_expression = det_ratio * (np.sum(s0) + self.a_alpha_prior) / (
                (len(self.s) - np.sum(s0) + self.a_beta_prior -np.float32(1)))
            exp = np.exp(exp_subtraction)
            ratio = pre_exp_expression * exp
            s0_prob = np.float32(1) / (np.float32(1) + ratio.item())
            s1_prob = np.float32(1) - s0_prob

        elif exp_subtraction > self.exp_subtraction_threshold:
            s1_prob = np.float32(1)
        elif exp_subtraction < -self.exp_subtraction_threshold:
            s1_prob = np.float32(0)

        changed = np.float32(0)

        from_0_to_1 = np.float32(0)
        from_1_to_0 = np.float32(0)

        if inserted_0:
            # put the new feature back into natural feature order
            s1 = self.s.copy()
            s1[int(s_index_to_sample)] = 1
            dest = int(np.sum(s1[:int(s_index_to_sample)]))
            xtx_1, xty_1 = move_front_back(xtx_1, xty_1, dest)

        if np.random.random() < s1_prob:
            if self.s[int(s_index_to_sample)] == 0:
                changed = np.float32(1)
                from_0_to_1 = np.float32(1)
            self.s[int(s_index_to_sample)] = np.float32(1)
        else:
            if self.s[int(s_index_to_sample)] == 1:
                changed = np.float32(1)
                from_1_to_0 = np.float32(1)
            self.s[int(s_index_to_sample)] = np.float32(0)

        return changed, from_0_to_1, from_1_to_0, xtx_0, xty_0, xtx_1, xty_1

    def sample_a(self, current_size, s, a_alpha_prior, a_beta_prior):
        a_alpha = current_size + a_alpha_prior
        a_beta = s.size - current_size + a_beta_prior
        self.a = np.random.beta(np.float32(a_alpha), np.float32(a_beta), size=1)[0]


    def sample_tau(self, current_size, w, tau_alpha_prior, tau_beta_prior):
        alpha = tau_alpha_prior + current_size / np.float32(2)
        w_contig = np.ascontiguousarray(w)
        beta = np.float32(1) / (np.float32(0.5) * (w_contig @ w_contig) + tau_beta_prior)
        gamma_draw = np.random.gamma(alpha, scale=beta, size=1).astype(np.float32)[0]
        self.tau = gamma_draw ** np.float32(-0.5)

    def sample_sigma(self, current_size, xs, w, y, sigma_alpha_prior, sigma_beta_prior):
        if current_size == 0:
            e = np.ascontiguousarray(y)
        else:
            w_contig = np.ascontiguousarray(w)
            xw = np.ascontiguousarray(xs @ w_contig)
            e = np.ascontiguousarray(y - xw)

        sigma_alpha = sigma_alpha_prior + xs.shape[0] / np.float32(2.0)
        e_dot = np.dot(e, e)
        sigma_beta = np.float32(1.0) / (np.float32(0.5) * e_dot + sigma_beta_prior)

        gamma_draw = np.random.gamma(sigma_alpha, scale=sigma_beta, size=1).astype(np.float32)[0]
        self.sigma = gamma_draw ** np.float32(-0.5)

    def sample_w(self, xs, y, tau, sigma, s_size, xtx, xty):
        k = int(s_size)

        if k == 0:
            self.w = np.zeros(1, dtype=np.float32)
            return np.zeros(1, dtype=np.float32)

        Q = (np.float32(1.0) / (tau * tau)) * np.eye(k) + (np.float32(1.0) / (sigma * sigma)) * xtx

        L = np.linalg.cholesky(Q).astype(np.float32)
        rhs = (np.float32(1.0) / (sigma * sigma)) * xty
        u = np.linalg.solve(L, rhs.astype(np.float32))
        mu = np.linalg.solve(L.T.astype(np.float32), u.astype(np.float32))
        z = np.random.normal(0.0, 1.0, k).astype(np.float32)

        # PATCH 2: Q = L L^T is the precision, so the noise with covariance Q^-1 is L^-T z
        v = np.linalg.solve(L.T.astype(np.float32), z.astype(np.float32))
        result = mu + v

        self.w = result
        return result

    def cycle(self, save=True):
        index_order = np.random.permutation(self.x.shape[1]).astype(np.float32)  # Numba-compatible
        xtx, xty = build_xtx_xty_pair(self.x, self.y, self.s)

        dummy = np.zeros((0, 0), dtype=np.float32)
        dummy_vec = np.zeros(0, dtype=np.float32)

        for i in range(index_order.shape[0]):
            s_index_to_sample = index_order[i]

            s_for_feature = self.s[int(s_index_to_sample)]

            if s_for_feature == 0:
                changed, from_0_to_1, from_1_to_0, xtx_0, xty_0, xtx_1, xty_1 = (
                    self.sample_s(s_index_to_sample, xtx_0=xtx, xty_0=xty, xtx_1=dummy, xty_1=dummy_vec, inserted_0=True))
            else:
                changed, from_0_to_1, from_1_to_0, xtx_0, xty_0, xtx_1, xty_1 = (
                    self.sample_s(s_index_to_sample, xtx_1=xtx, xty_1=xty, xtx_0=dummy, xty_0=dummy_vec, inserted_0=False))
            if changed:
                if from_0_to_1:
                    xtx = xtx_1
                    xty = xty_1
                else:
                    xtx = xtx_0
                    xty = xty_0

        # PATCH 1: w, a, tau and sigma are redrawn once in every sweep, whether or not a switch flipped
        self.sample_all_besides_s(xtx, xty)

        if save:
            self.a_chain[self.t] = self.a
            self.s_chain[:, self.t] = self.s
            self.tau_chain[self.t] = self.tau
            self.sigma_chain[self.t] = self.sigma
            self.w_chain[:, self.t] = insert_masked(self.s, self.w)
            self.w_dim_chain[self.t] = np.sum(self.s, dtype=np.int32)
            self.t += 1

    def sample_all_besides_s(self, xtx, xty):
        xs = get_xs(self.s, self.x)
        current_size = np.float32(np.sum(self.s))
        self.sample_w(xs, self.y, self.tau, self.sigma, current_size, xtx, xty)
        self.sample_a(current_size, self.s, self.a_alpha_prior, self.a_beta_prior)
        self.sample_tau(current_size, self.w, self.tau_alpha_prior, self.tau_beta_prior)
        self.sample_sigma(current_size, xs, self.w, self.y, self.sigma_alpha_prior, self.sigma_beta_prior)


    def run_cycles(self, n_steps, start_collect, sample_geweke=False):
        for i in range(n_steps):
            if sample_geweke:
                self.sample_geweke()

            self.cycle(True)

            # collect internal predictions only in the last 10 %
            if self.t >= start_collect:
                self.compute_internal_y_prediction()
