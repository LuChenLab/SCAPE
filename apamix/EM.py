import sys
import random
import traceback
import pandas as pd
import numpy as np


from loguru import logger
from scipy.stats import norm, lognorm
from rpy2.robjects.packages import STAP
from rpy2.robjects.numpy2ri import numpy2rpy

from apamix.mix_utils import *
from utils.utils import dotdict


class EM:
    @staticmethod
    def lik_f0(L, max_LA, log=True):

        pl_s = px = 1 / L
        pr_s = 1 / max_LA
        res = px * pl_s * pr_s
        if log:
            return np.log(res)
        else:
            return res

    @staticmethod
    def lik_r_s(r_arr,
                s,
                log=False):
        res = (r_arr <= s) / s
        if log:
            return np.log(res)
        else:
            return res

    @staticmethod
    def lik_x_st(x_arr,
                 s,
                 theta,
                 mu_f,
                 sigma_f,
                 log=False):
        """
        calculate p(l|x, theta)
        """
        if log:
            res = lognorm.pdf(x_arr,
                              theta + s + 1 - mu_f,
                              sigma_f)
        else:
            res = norm.pdf(x_arr,
                           theta + s + 1 - mu_f,
                           sigma_f)
        return res

    @staticmethod
    def lik_l_xt(x_arr,
                 l_arr,
                 theta,
                 log=False):
        # calculate the UTR part length by given the pa site
        utr_len_arr = theta - x_arr + 1
        # neg value can't support the given theta
        valid_inds = l_arr <= utr_len_arr

        res = valid_inds + 0

        if len(valid_inds) == 1:
            # for single Flase array
            if all(valid_inds):
                res = np.array(1/utr_len_arr[valid_inds])
        else:
            res[valid_inds] = np.divide(1, utr_len_arr[valid_inds])

        if any(np.isinf(res)):
            ValueError("res contains Inf in lik_l_xt calculation")

        if log:
            return np.log(res)
        else:
            return res

    @staticmethod
    def lik_lsr_t(x_arr,
                  l_arr,
                  r_arr,
                  s_arr,
                  theta,
                  mu_f,
                  sigma_f,
                  pmf_s_dis_arr,
                  s_dis_arr,
                  log=False):
        s_len = len(s_dis_arr)
        res = np.zeros(len(x_arr))

        oth_inds = np.isnan(s_arr)
        valid_inds = np.invert(oth_inds)

        n_valid_inds = np.sum(valid_inds)
        n_other_inds = np.sum(oth_inds)

        if n_valid_inds > 0:
            tmp1 = EM.lik_x_st(
                x_arr[valid_inds],
                s_arr[valid_inds],
                theta,
                mu_f,
                sigma_f,
                log=log
            )

            tmp2 = EM.lik_l_xt(
                x_arr[valid_inds],
                l_arr[valid_inds],
                theta,
                log=log
            )

            if log:
                res[valid_inds] = tmp1 + tmp2
            else:
                res[valid_inds] = tmp1 * tmp2

        if not any(oth_inds):
            return (res)

        tmp_sum = np.zeros(n_other_inds)

        for i in np.arange(s_len):
            tmp1 = EM.lik_r_s(r_arr[oth_inds], s_dis_arr[i])
            tmp2 = EM.lik_x_st(x_arr[oth_inds],
                               s_dis_arr[i],
                               theta,
                               mu_f,
                               sigma_f,
                               log=log)

            tmp3 = EM.lik_l_xt(x_arr[oth_inds],
                               l_arr[oth_inds],
                               theta,
                               log=log)
            tmp_sum += tmp1 * tmp2 * tmp3 * pmf_s_dis_arr[i]

        if log:
            res[oth_inds] = np.log(tmp_sum)
        else:
            res[oth_inds] = tmp_sum
        return res


    @staticmethod

    def lik_lsr_t0(x_arr,
                   l_arr,
                   r_arr,
                   s_arr,
                   theta,
                   pmf_s_dis_arr,
                   s_dis_arr,
                   log=False):
        """
        calculate p(x,l,r|theta) or log(p(x,l,r,|theta)),
        theta is known
        """

        s_len = len(s_dis_arr)
        try:
            res = np.zeros(len(x_arr))
            oth_inds = np.isnan(s_arr)
            valid_inds = np.invert(oth_inds)
        except TypeError:
            x_arr = np.array([x_arr])
            l_arr = np.array([l_arr])
            r_arr = np.array([r_arr])
            res = np.zeros(1)
            oth_inds = np.array([np.isnan(s_arr)])
            valid_inds = np.array(np.invert(oth_inds))

        n_valid_inds = np.sum(valid_inds)
        n_other_inds = np.sum(oth_inds)

        if n_valid_inds > 0:
            tmp1 = EM.lik_l_xt(x_arr[valid_inds],
                               l_arr[valid_inds],
                               theta,
                               log=log)

            if log:
                tmp2 = 0
                res[valid_inds] = tmp1 + tmp2
            else:
                tmp2 = 1
                res[valid_inds] = tmp1 * tmp2

        if not any(oth_inds):
            return res

        tmp_sum = np.zeros(n_other_inds)

        for i in np.arange(s_len):
            tmp1 = EM.lik_r_s(r_arr[oth_inds], s_dis_arr[i])
            tmp2 = 1
            tmp3 = EM.lik_l_xt(x_arr[oth_inds],
                               l_arr[oth_inds],
                               theta,
                               log=False)
            tmp_sum += tmp1 * tmp2 * tmp3 * pmf_s_dis_arr[i]

        if log:
            res[oth_inds] = np.log(tmp_sum)
        else:
            res[oth_inds] = tmp_sum

        return res

    @staticmethod
    def em_algo(ws,
                theta_win_mat,
                n_frag,
                pre_alpha_arr,
                pre_beta_arr,
                all_theta,
                init_beta_arr,
                all_win_log_lik_mat_list,
                pre_init_alpha_arr,
                min_pa_gap,
                unif_log_lik,
                min_theta,
                L,
                nround=50,
                verbose=False):
        lb = -np.Inf
        lb_arr = np.repeat(np.NaN, nround)
        K = len(ws) - 1

        alpha_arr = theta_win_mat.res_alpha_arr
        beta_arr = theta_win_mat.res_beta_arr

        if pre_alpha_arr:
            pre_flag_arr = np.isin(alpha_arr, pre_init_alpha_arr)

        if pre_alpha_arr and pre_beta_arr:
            tmp_pre_init_alpha_arr = pre_init_alpha_arr[
                np.invert(
                    np.isnan(pre_beta_arr)
                )
            ]
            tmp_arr = gen_k_arr(
                len(alpha_arr) - len(tmp_pre_init_alpha_arr), nround
            )
            opt_inds = np.arange(len(alpha_arr))[
                np.isin(alpha_arr,
                        tmp_pre_init_alpha_arr,
                        invert=True)
            ]
            k_arr = tmp_arr.copy()
            for i in range(len(opt_inds)):
                k_arr[tmp_arr == i] = opt_inds[i]

        elif pre_alpha_arr and not pre_beta_arr:
            k_arr = gen_k_arr(len(alpha_arr), nround)

        elif not pre_alpha_arr:
            pre_init_alpha_arr = None
            k_arr = gen_k_arr(len(alpha_arr), nround)
        log_zmat = np.zeros(
            (n_frag, K + 1)
        )

        for k in range(K + 1):
            log_zmat = EM.cal_z_k(
                all_win_log_lik_mat_list,
                alpha_arr,
                beta_arr,
                ws,
                k,
                log_zmat,
                all_theta,
                init_beta_arr,
                unif_log_lik
            )

        for i in range(nround):
            if verbose:
                logger.debug(f'iteration={i+1}, lb={lb}')

            # e step
            log_zmat = EM.cal_z_k(
                all_win_log_lik_mat_list,
                alpha_arr,
                beta_arr,
                ws,
                k_arr[i],
                log_zmat,
                all_theta,
                init_beta_arr,
                unif_log_lik
            )

            Z = EM.norm_z(log_zmat)

            res = EM.mstep(
                alpha_arr,
                beta_arr,
                ws,
                Z,
                k_arr[i],
                all_theta,
                min_theta,
                L,
                init_beta_arr,
                all_win_log_lik_mat_list,
                pre_init_alpha_arr=pre_init_alpha_arr
            )

            alpha_arr = res.alpha_arr
            beta_arr = res.beta_arr
            ws = res.ws

            lb_new = EM.elbo(log_zmat, Z)
            lb_arr[i] = lb_new
            if lb_new == -np.Inf:
                lb = -np.Inf
                break

            if abs(lb_new - lb) < abs(1e-6 * lb):
                break
            else:
                lb = lb_new
        if i == nround and verbose:
            logger.debug(f'Run all {i} iterations. lb={lb}')
        elif verbose:
            logger.debug(f'Run all {i} iterations. lb={lb}')
        else:
            pass

        bic = EM.cal_bic(log_zmat, Z)

        if verbose:
            logger.debug(
                f'''
                bic={bic}\n
                estimated ws: {ws} \n
                estimated alpha: {list(alpha_arr)} \n
                estimated beta: {list(beta_arr)} \n
                ''')

        if pre_alpha_arr:
            oth_inds = np.where(np.invert(pre_flag_arr))
            for i in oth_inds:
                if any(abs(alpha_arr[i] - pre_alpha_arr) < min_pa_gap):
                    logger.warning(f'alpha: {alpha_arr[i]} is whith {min_pa_gap} distance')
                    return None

        lb_arr = lb_arr[np.invert(np.isnan(lb_arr))]
        res = dotdict({
            'ws': ws,
            'alpha_arr': alpha_arr,
            'beta_arr': beta_arr,
            'lb_arr': lb_arr,
            'bic': bic
        })
        return res

    @staticmethod
    def mstep(alpha_arr,
              beta_arr,
              ws,
              Z,
              k,
              all_theta,
              min_theta,
              L,
              init_beta_arr,
              all_win_log_lik_mat_list,
              pre_init_alpha_arr=None):
        K = len(ws) - 1
        new_ws = EM.maximize_ws(Z)

        if np.isin(alpha_arr[k], pre_init_alpha_arr):
            valid_inds = np.where(all_theta == alpha_arr[k])[0]
        else:
            tmp_min_theta = min_theta if k == 0 else alpha_arr[k - 1] + beta_arr[k - 1]
            tmp_max_theta = L if k == K - 1 else alpha_arr[k + 1]
            valid_inds = np.where(np.logical_and(
                all_theta >= tmp_min_theta,
                all_theta <= tmp_max_theta
            ))[0]
        n_winsize = len(init_beta_arr)
        res_list = []
        loglik_arr = np.repeat(-np.Inf, n_winsize)

        for i in range(n_winsize):
            res_list.append(
                EM.maximize_win_k(valid_inds, i, Z[:, k], all_theta, init_beta_arr, all_win_log_lik_mat_list)
            )
            loglik_arr[i] = res_list[i].loglik
        max_ind = np.argmax(loglik_arr)
        res = res_list[max_ind]

        alpha_arr[k] = res.alpha
        beta_arr[k] = res.beta
        res = dotdict({
            'alpha_arr': alpha_arr,
            'beta_arr': beta_arr,
            'ws': new_ws
        })
        return res

    @staticmethod
    def maximize_win_k(theta_inds, beta_ind, zk, all_theta, init_beta_arr, all_win_log_lik_mat_list):
        n_input_theta = len(theta_inds)
        valid_inds = zk > 0

        beta = init_beta_arr[beta_ind]
        all_win_log_lik_mat = all_win_log_lik_mat_list[beta_ind]

        tmp_lik_arr = np.zeros(n_input_theta)
        for i in range(n_input_theta):
            tmparr = all_win_log_lik_mat[valid_inds, theta_inds[i]]
            tmp_lik_arr[i] = np.sum(zk[valid_inds] * tmparr)
        max_ind = np.argmax(tmp_lik_arr)
        max_loglik = tmp_lik_arr[max_ind]
        alpha = all_theta[theta_inds[max_ind]]
        res = dotdict({
            "alpha": alpha,
            "beta": beta,
            "loglik": max_loglik
        })
        return res

    @staticmethod
    def maximize_ws(Z):
        ws = Z.sum(axis=0) / Z.shape[0]
        return ws

    @staticmethod
    def exp_log_lik(log_zmat, Z):
        ZZ = Z * log_zmat
        ZZ[Z == 0] = 0
        return np.sum(ZZ)

    @staticmethod
    def elbo(log_zmat, Z):
        LZ = Z.copy()
        LZ[Z != 0] = np.log(Z[Z != 0])
        entropy = -1 * Z * LZ
        lb = EM.exp_log_lik(log_zmat, Z) + np.sum(entropy)
        if np.isnan(lb):
            raise ValueError('lower bounder is na.')
        return lb

    @staticmethod
    def cal_z_k(all_win_log_lik_mat_list,
                alpha_arr,
                beta_arr,
                ws,
                k,
                log_zmat,
                all_theta,
                init_beta_arr,
                unif_log_lik):
        K = len(ws) - 1
        tmp_len = len(all_win_log_lik_mat_list)

        if k < K:
            tmp_ind = np.arange(tmp_len)[np.array(init_beta_arr == beta_arr[k])]
            all_win_log_lik_mat = all_win_log_lik_mat_list[tmp_ind[0]]
            log_zmat[:, k] = np.log(ws[k]) + \
                             all_win_log_lik_mat[:, all_theta == alpha_arr[k]].flat
        else:
            log_zmat[:, k] = np.log(ws[k]) + unif_log_lik
        return log_zmat

    @staticmethod
    def fixed_inference(alpha_arr,
                        beta_arr,
                        n_frag,
                        nround,
                        all_theta,
                        all_win_log_lik_mat_list,
                        init_beta_arr,
                        unif_log_lik,
                        verbose=False):
        lb = -np.Inf
        lb_arr = np.repeat(np.NaN, nround)
        K = len(alpha_arr)

        ws = np.random.uniform(size=K + 1) + 1
        ws[K] += -1
        ws = ws / ws.sum()

        k_arr = gen_k_arr(K, nround)
        log_zmat = np.zeros(
            (n_frag, K + 1)
        )

        for k in range(K + 1):
            log_zmat = EM.cal_z_k(
                all_win_log_lik_mat_list,
                alpha_arr,
                beta_arr,
                ws,
                k,
                log_zmat,
                all_theta,
                init_beta_arr,
                unif_log_lik
            )

        for i in range(nround):
            if verbose:
                logger.debug(f'iteration={i}, lb={lb}')

            # e step
            log_zmat = EM.cal_z_k(
                all_win_log_lik_mat_list,
                alpha_arr,
                beta_arr,
                ws,
                k_arr[i],
                log_zmat,
                all_theta,
                init_beta_arr,
                unif_log_lik
            )

            Z = EM.norm_z(log_zmat)

            ws = EM.maximize_ws(Z)

            lb_new = EM.elbo(log_zmat, Z)
            lb_arr[i] = lb_new
            if lb_new == -np.Inf:
                lb = -np.Inf
                break

            if abs(lb_new - lb) < abs(1e-6 * lb):
                break
            else:
                lb = lb_new
        if i == nround and verbose:
            logger.debug(f'Run all {i} iterations. lb={lb}')
        elif verbose:
            logger.debug(f'Run all {i} iterations. lb={lb}')
        else:
            pass

        bic = EM.cal_bic(log_zmat, Z)

        label = log_zmat.argmax(axis=1)
        if verbose:
            logger.debug(
                f'''
                bic={bic}\n
                estimated ws: {ws} \n
                estimated alpha: {list(alpha_arr)} \n
                estimated beta: {list(beta_arr)} \n
                ''')

        lb_arr = lb_arr[np.invert(np.isnan(lb_arr))]
        res = dotdict({
            'ws': ws,
            'alpha_arr': alpha_arr,
            'beta_arr': beta_arr,
            'lb_arr': lb_arr,
            'bic': bic,
            'label': label
        })
        return res

    @staticmethod
    def norm_z(Z):
        Z = np.exp(Z - np.max(Z, axis=1)[:, None])
        Z = Z / Z.sum(axis=1)[:, None]
        return Z

    @staticmethod
    def cal_bic(log_zmat, Z):
        N = np.shape(Z)[0]
        K = np.shape(Z)[1] - 1

        # the smaller bic, the better model
        res = -2 * EM.exp_log_lik(log_zmat, Z) + (3 * K + 1) * np.log(N)
        return res
