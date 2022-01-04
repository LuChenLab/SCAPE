import sys
import random
import traceback
import pickle

import pandas as pd
import numpy as np

from loguru import logger
from scipy.stats import norm, lognorm
from rpy2.robjects.packages import STAP
from rpy2.robjects.numpy2ri import numpy2rpy

from apamix.EM import EM
from apamix.mix_utils import *
from utils.utils import dotdict, dict_list


class APA:

    def __init__(self,
                 n_max_apa,
                 n_min_apa,
                 r1_utr_st_arr,
                 r1_len_arr,
                 r2_len_arr,
                 polya_len_arr,
                 pa_site_arr,
                 utr_len,
                 LA_dis_arr=None #  [10, 30, 50, 70, 90, 110, 130],
                 pmf_LA_dis_arr=None # [309912, 4107929, 802856, 518229, 188316, 263208, 101],
                 min_LA=20,
                 max_LA=150,
                 mu_f=300,
                 sigma_f=50,
                 min_ws=0.01,
                 min_pa_gap=100,
                 max_beta=70,
                 pre_alpha_arr=None,
                 pre_beta_arr=None,
                 theta_step=9,
                 cb=None,
                 umi=None,
                 fixed_inference_flag=False,
                 pre_init_theta_win_mat=None,
                 region_name=None,
                 verbose=False
                 ):

        self.n_max_apa = n_max_apa
        self.n_min_apa = n_min_apa
        self.r1_utr_st_arr = r1_utr_st_arr
        self.r1_len_arr = r1_len_arr
        self.min_theta = min(self.r1_len_arr)
        self.r2_len_arr = r2_len_arr
        self.verbose = verbose

        self.polya_len_arr = np.array(polya_len_arr)
        self.pa_site_arr = np.array(pa_site_arr)
        self.L = utr_len
        self.min_LA = min_LA
        self.max_LA = max_LA

        self.mu_f = mu_f
        self.sigma_f = sigma_f
        self.min_ws = min_ws

        self.theta_step = theta_step
        self.n_frag = len(self.r1_len_arr)
        self.min_pa_gap = min_pa_gap

        self.region_name = region_name

        if LA_dis_arr:
            self.s_dis_arr = np.array(LA_dis_arr)
            self.pmf_s_dis_arr = np.array(pmf_LA_dis_arr)
        else:
            self.s_dis_arr = np.arange(self.min_LA, self.max_LA, 10)
            self.pmf_LA_dis_arr = np.repeat(1 / len(self.s_dis_arr), len(self.s_dis_arr))
        self.pmf_s_dis_arr = self.pmf_s_dis_arr / sum(self.pmf_s_dis_arr)
        self.nround = 50

        # For dropout single cell rna seq
        self.cb = cb
        self.umi = umi

        self.pre_alpha_arr = np.array(pre_alpha_arr) if pre_alpha_arr else pre_alpha_arr
        self.pre_beta_arr = np.array(pre_beta_arr) if pre_beta_arr else pre_beta_arr
        self.max_beta = max_beta
        self.all_theta = np.arange(self.min_theta, self.L, self.theta_step)
        self.n_all_theta = len(self.all_theta)

        self.__check()  # Check the alpha and beta arr information

        self.fixed_inference_flag = fixed_inference_flag
        self.pre_init_theta_win_mat = pre_init_theta_win_mat

        self.unif_log_lik = EM.lik_f0(self.L, self.max_LA, log=True)
        self.__init_run()

        self.data_wins = self.__split_data()


    def __check(self):
        if all(np.isnan(self.pa_site_arr)):
            self.pa_site_arr = np.repeat(np.NaN, self.n_frag)

        if all(np.isnan(self.polya_len_arr)):
            self.polya_len_arr = np.repeat(np.NaN, self.n_frag)

        # To check the alpha and beta arr was legal
        n_pre_pa_sites = 0
        if self.pre_alpha_arr:
            self.min_theta = min(min(self.r1_len_arr),
                                 min(self.pre_alpha_arr))
            self.pre_alpha_arr.sort()
            if max(self.pre_alpha_arr) > self.L:
                raise ValueError(
                    f'''
                    The max of pre_alpha_arr was out of UTR length.
                    '''
                )
            n_pre_pa_sites = len(self.pre_alpha_arr)

            self.pre_init_alpha_arr = replace_with_closest(self.all_theta, self.pre_alpha_arr)
        else:
            self.pre_init_alpha_arr = None
        if self.pre_beta_arr:
            if n_pre_pa_sites != len(self.pre_beta_arr):
                raise ValueError(
                    f'''
                    The num of pre_pa_sites wasn't same with pre_beta_arr.
                    '''
                )
            if not all(self.pre_beta_arr):
                raise ValueError(
                    f'''
                    All pre_beta_arr was NaN.
                    '''
                )
            if np.max(self.pre_beta_arr) > self.max_beta:
                raise ValueError(
                    f'''
                    The max of pre_beta_arr was more than {self.max_beta}
                    '''
                )

        if n_pre_pa_sites > 0:
            if n_pre_pa_sites > self.n_max_apa:
                logger.warning(
                    f'''
                    n_max_apa: {self.n_max_apa} change to n_pre_pa_sites: {n_pre_pa_sites}
                    '''
                )
                self.n_max_apa = n_pre_pa_sites

            if n_pre_pa_sites > self.n_min_apa:
                logger.warning(
                    f'''
                    n_min_apa: {self.n_max_apa} change to n_pre_pa_sites: {n_pre_pa_sites}
                    '''
                )
                self.n_min_apa = n_pre_pa_sites

        if not self.max_beta >= self.theta_step:
            raise ValueError(
                f'''
                max_beta: {self.max_beta} wasn't more than {self.theta_step}
                '''
            )
        if self.n_frag != len(self.r1_utr_st_arr) or self.n_frag != len(self.r2_len_arr):
            raise ValueError(
                f'''
                n_frag ({self.n_frag}) wasn't same to r1_utr_st_arr ({self.r1_utr_st_arr})
                '''
            )

        if self.L < self.min_theta:
            raise ValueError(
                f'''
                L ({self.L}) wasn't more than min_theta ({self.min_theta})
                '''
            )

        return self

    def __init_run(self):
        self.lb_arr = np.repeat(-np.Inf, self.n_max_apa)
        self.bic_arr = np.repeat(np.Inf, self.n_max_apa)
        self.ws_flag_arr = np.repeat(False, self.n_max_apa)
        self.res_lst = []

        self.all_theta_lik_mat = np.zeros((self.n_frag, self.n_all_theta), dtype='float64')

        oth_inds = np.isnan(self.pa_site_arr)
        pa_inds = np.invert(oth_inds)

        for i in np.arange(self.n_all_theta):
            self.all_theta_lik_mat[oth_inds, i] = EM.lik_lsr_t(
                self.r1_utr_st_arr[oth_inds],
                self.r1_len_arr[oth_inds],
                self.r2_len_arr[oth_inds],
                self.polya_len_arr[oth_inds],
                self.all_theta[i],
                self.mu_f,
                self.sigma_f,
                self.pmf_s_dis_arr,
                self.s_dis_arr
            )

        if any(pa_inds):
            for i in np.where(pa_inds)[0]:
                tmp_ind = np.argmin(abs(self.all_theta - self.pa_site_arr[i]))
                self.all_theta_lik_mat[i, tmp_ind] = EM.lik_lsr_t0(
                    self.r1_utr_st_arr[i],
                    self.r1_len_arr[i],
                    self.r2_len_arr[i],
                    self.polya_len_arr[i],
                    self.all_theta[tmp_ind],
                    self.pmf_s_dis_arr,
                    self.s_dis_arr
                )

        self.init_beta_arr = np.arange(
            self.theta_step,
            np.ceil(
                self.max_beta / self.theta_step
            ) * self.theta_step + 1.0,
            self.theta_step)

        self.all_win_log_lik_mat_list = []
        for i in range(len(self.init_beta_arr)):
            self.all_win_log_lik_mat_list.append(
                self.kernel_smooth(
                    self.init_beta_arr[i]
                )
            )

        return self

    def __split_data(self):
        n_frag = len(self.r1_utr_st_arr)
        coverage_cnt = np.zeros(self.L)

        for i in range(n_frag):
            tmp_inds = self.r1_utr_st_arr[i] - 1 + np.arange(self.r1_len_arr[i])
            coverage_cnt[tmp_inds] += 1


        ksmooth = '''run = function(L, coverage_cnt, theta_step){
                ks_res = stats::ksmooth(
                        seq(L),
                        coverage_cnt,
                        kernel = 'normal',
                        bandwidth = 5 * theta_step,
                        x.points = seq(L)
                    )
                return(ks_res$y)
            }
            '''
        ksmooth = STAP(ksmooth, 'run')

        ks_res = np.array(
            ksmooth.run(
                self.L,
                numpy2rpy(coverage_cnt),
                self.theta_step
            )
        )
        ks_res_x = np.arange(self.L) + 1
        sign_arr = np.sign(np.diff(np.concatenate(([-1], ks_res))))

        if sign_arr[0] != 1:
            raise ValueError('First sign value not 1, split_data func')

        if not any(sign_arr == 0):
            pass
        else:
            st_arr, lengths, tmp_vals = rle(sign_arr)
            en_arr = np.cumsum(lengths)
            tmp_n = len(tmp_vals)
            zero_inds = np.where(tmp_vals == 0)

            for _, i in np.ndenumerate(zero_inds):
                st = st_arr[i]
                en = en_arr[i]
                if i == tmp_n - 1 or en - st <= 1:
                    sign_arr[st:en] = tmp_vals[i - 1]
                else:
                    mid = np.int(np.round((st + en + 1) / 2))
                    sign_arr[st:mid] = tmp_vals[i - 1]
                    sign_arr[mid:en] = tmp_vals[i + 1]
        st_arr, lengths, tmp_vals = rle(sign_arr)
        chng = np.cumsum(lengths)
        mode_arr = ks_res_x[chng[np.arange(len(chng), step=2)] - 1]
        n_mode = len(mode_arr)
        boarder_arr = ks_res_x[chng[np.arange(1, len(chng), step=2)] - 1][:n_mode - 1]
        st_arr = np.concatenate(([1], boarder_arr + 1))
        en_arr = np.concatenate((boarder_arr, [self.L]))
        ws_arr = np.zeros(n_mode) + 1

        for i in range(n_mode):
            ws_arr[i] = sum(coverage_cnt[st_arr[i]:en_arr[i]])
        ws_arr = ws_arr / ws_arr.sum()

        return pd.DataFrame({
            "st_arr": st_arr,
            "en_arr": en_arr,
            "ws_arr": ws_arr,
            "mode_arr": mode_arr
        }
        )

    def kernel_smooth(self, beta):
        n_step = beta / self.theta_step
        n_bd_step = int(np.ceil(n_step * 3))
        weights = norm.pdf(np.arange(-n_bd_step, n_bd_step), 0, n_step)
        weights = weights / np.sum(weights)

        n_all_theta = len(self.all_theta)

        if n_all_theta != self.all_theta_lik_mat.shape[1]:
            raise ValueError(
                f'''
                n_all_theta: {n_all_theta} wasn't equal to column of all_theta_lik_mat
                '''
                )

        all_win_log_lik_mat = np.zeros((self.n_frag,
                                        n_all_theta))
        for i in range(n_bd_step):
            tmp_inds = np.arange(np.int(i - n_bd_step), np.int(i + n_bd_step))

            tmp_inds = tmp_inds[tmp_inds >= 0]

            tmp_weight = weights[tmp_inds]
            tmp_weight = tmp_weight / np.sum(tmp_weight)
            tmp_weight = np.reshape(tmp_weight, (-1, len(tmp_weight)))
            all_win_log_lik_mat[:, i] = np.log(
                np.sum(
                    tmp_weight * self.all_theta_lik_mat[:, tmp_inds],
                    axis=1
                )
            )
        for i in range(n_bd_step, n_all_theta - n_bd_step):
            tmp_inds = np.arange(
                np.int(i - n_bd_step),
                np.int(i + n_bd_step)
            )
            tmp_weight = weights
            all_win_log_lik_mat[:, i] = np.log(
                np.sum(
                    tmp_weight * self.all_theta_lik_mat[:, tmp_inds],
                    axis=1
                )
            )
        for i in range(n_all_theta - n_bd_step, n_all_theta):
            tmp_inds = np.arange(
                np.int(i - n_bd_step),
                np.int(i + n_bd_step)
            )
            tmp_weight = weights[tmp_inds < n_all_theta]
            tmp_weight = tmp_weight / np.sum(tmp_weight)
            tmp_weight = np.reshape(tmp_weight, (-1, len(tmp_weight)))

            tmp_inds = tmp_inds[tmp_inds < n_all_theta]

            all_win_log_lik_mat[:, i] = np.log(
                np.sum(
                    tmp_weight * self.all_theta_lik_mat[:, tmp_inds],
                    axis=1
                )

            )

        return all_win_log_lik_mat

    def init_theta_win(self,
                       n_apa):
        def sample_theta(k_arr, all_theta, data_wins, mode='full'):
            res_arr = np.repeat(0.0, len(k_arr))
            for i in range(len(k_arr)):
                iw = k_arr[i]
                if mode == 'full':
                    tmp_theta_arr = \
                        all_theta[
                            np.logical_and(all_theta >= data_wins.st_arr[iw],
                                           all_theta < data_wins.en_arr[iw])
                        ]
                elif mode == 'left':
                    tmp_theta_arr = \
                        all_theta[
                            np.logical_and(all_theta >= data_wins.st_arr[iw],
                                           all_theta < data_wins.mode_arr[iw])
                        ]
                elif mode == 'right':
                    tmp_theta_arr = \
                        all_theta[
                            np.logical_and(
                                all_theta >= data_wins.mode_arr[iw],
                                all_theta < data_wins.en_arr[iw]
                            )
                        ]
                else:
                    raise ValueError(f'Unknown mode: {mode}')

                if len(tmp_theta_arr) > 0:
                    res_arr[i] = np.random.choice(
                        a=tmp_theta_arr,
                        size=1,
                        replace=False
                    )
                else:
                    res_arr[i] = -1
            return np.sort(res_arr)

        # special treatment for the second sampling
        def sample_theta1(k_arr, n_data_win, all_theta, data_wins):
            tmp_ind = (k_arr == (n_data_win))
            if tmp_ind.sum() > 0:
                theta_arr1 = sample_theta([n_data_win-1], all_theta, data_wins, mode='full')
            else:
                theta_arr1 = np.array([])
            tmp_ind = (k_arr != (n_data_win))
            if tmp_ind.sum() > 0:
                theta_arr2 = sample_theta(k_arr[tmp_ind], all_theta, data_wins, mode='left')
            else:
                theta_arr2 = np.array([])
            return np.append(theta_arr1, theta_arr2)

        def get_data_win_lab(alpha_arr, data_wins):
            n_data_win = len(data_wins.st_arr)
            left = 0
            right = data_wins.mode_arr[0]
            lab_arr = alpha_arr.copy()
            for i in range(n_data_win):
                lab_arr[np.logical_and(alpha_arr >= left, alpha_arr < right)] = i
                left = data_wins.mode_arr[i]
                right = data_wins.mode_arr[i + 1] if i < n_data_win - 1 else None

            lab_arr[alpha_arr >= left] = n_data_win
            return lab_arr

        # process pre_alpha_arr and pre_beta_arr
        pre_para_flag = True if self.pre_alpha_arr else False
        if pre_para_flag:
            pre_init_alpha_arr = replace_with_closest(
                self.all_theta,
                self.pre_alpha_arr
            )
            pre_init_beta_arr = np.random.choice(
                self.init_beta_arr[np.arange(np.int(len(self.init_beta_arr) / 2 + 1))],
                size=len(self.pre_alpha_arr),
                replace=True
            )
            if self.pre_beta_arr:
                tmp_inds = np.invert(np.isnan(self.pre_beta_arr))
                pre_init_beta_arr[tmp_inds] = replace_with_closest(
                    self.init_beta_arr,
                    self.pre_beta_arr[tmp_inds]
                )

        if pre_para_flag and len(self.pre_alpha_arr) == n_apa:
            return pd.DataFrame(
                {'res_alpha_arr': pre_init_alpha_arr,
                 'res_beta_arr': pre_init_beta_arr}
            )

        k_big_inds = np.where(self.data_wins.ws_arr > 0.1)[0]
        n_data_win = len(self.data_wins.ws_arr)
        k1 = len(k_big_inds)

        n_max_trial = 200
        for i in range(n_max_trial):

            # try sample function
            if n_data_win >= n_apa:

                k_arr = np.sort(
                    np.random.choice(
                        a=n_data_win,
                        size=n_apa,
                        p=self.data_wins.ws_arr,
                        replace=False
                    )
                )

                tmp_alpha_arr = sample_theta(k_arr, self.all_theta, self.data_wins, mode='right')
            elif n_data_win >= n_apa - k1:
                theta_arr1 = sample_theta(np.arange(n_data_win), self.all_theta, self.data_wins, mode='right')
                if k1 == 1:
                    k_arr = k_big_inds + 1
                else:
                    ws_tmp = self.data_wins.ws_arr[k_big_inds]
                    k_arr = np.random.choice(
                        a=k_big_inds + 1,
                        size=n_apa - n_data_win,
                        p=ws_tmp/ws_tmp.sum(),
                        replace=False
                    )
                theta_arr2 = sample_theta1(k_arr, n_data_win, self.all_theta, self.data_wins)
                tmp_alpha_arr = np.sort(
                    np.append(theta_arr1, theta_arr2)
                )
            else:
                theta_arr1 = sample_theta(
                    np.arange(n_data_win),
                    self.all_theta,
                    self.data_wins,
                    mode='right'
                )
                theta_arr2 = sample_theta1(k_big_inds + 1, n_data_win, self.all_theta, self.data_wins)

                k_arr = np.random.choice(
                    a=n_data_win,
                    size=n_apa - n_data_win - k1,
                    p=self.data_wins.ws_arr
                )

                theta_arr3 = sample_theta(k_arr, self.all_theta, self.data_wins, mode='full')

                tmp_alpha_arr = np.sort(
                    np.concatenate(
                        (theta_arr1, theta_arr2, theta_arr3)
                    )
                )

            flag1 = all(np.diff(tmp_alpha_arr) >= self.min_pa_gap)
            flag2 = all(tmp_alpha_arr > 0)

            if flag1 and flag2:
                break

        if i == n_max_trial -1:
            raise ValueError(f'Failed to generate valid theta after {n_max_trial} trial.')

        if pre_para_flag:
            sam_lab = get_data_win_lab(tmp_alpha_arr, self.data_wins)
            pre_lab = get_data_win_lab(self.pre_alpha_arr, self.data_wins)
            n_to_rm = 0
            for i in range(len(pre_lab)):
                smp_inds = np.logical_and(
                    sam_lab == pre_lab[i],
                    sam_lab > 0
                )
                smp_inds_flag = np.sum(smp_inds)
                if smp_inds_flag == 0:
                    n_to_rm += 1
                elif smp_inds_flag == 1:
                    sam_lab[smp_inds] = 0
                else:
                    tmp_ind = np.argmin(
                        abs(
                            tmp_alpha_arr[smp_inds] - pre_alpha_arr[i]
                        )
                    )
                    sam_lab[smp_inds[tmp_ind]] = 0
            if n_to_rm > 0:
                smp_inds = sam_lab > 0
                if np.sum(smp_inds) > 1:
                    raise ValueError('spm_inds out of range, continue next trial')

                tmp_inds = np.random.choice(
                    a=smp_inds,
                    size=n_to_rm
                )
                sam_lab[tmp_inds] = 0

            valid_inds = sam_lab > 0
            res_alpha_arr = np.sort(np.concatenate(
                (tmp_alpha_arr[valid_inds],
                 self.pre_init_alpha_arr)
            ))
            res_beta_arr = np.random.choice(
                a=self.init_beta_arr[np.arange(int(len(self.init_beta_arr) / 2 + 1))],
                size=n_apa
            )
        else:
            res_alpha_arr = tmp_alpha_arr
            res_beta_arr = np.random.choice(
                a=self.init_beta_arr[np.arange(int(len(self.init_beta_arr) / 2 + 1))],
                size=n_apa
            )

        return pd.DataFrame(
            {
                'res_alpha_arr': res_alpha_arr,
                'res_beta_arr': res_beta_arr
            }
        )

    def em_optim(self,
                 n_apa,
                 n_trial=5,
                 verbose=False
                 ):
        lb_arr = np.repeat(-np.Infinity, n_trial)
        bic_arr = np.repeat(np.Infinity, n_trial)
        res_list = [None for i in range(n_trial)]

        for i in range(n_trial):
            try:
                # last component is for uniform component
                ws = np.random.uniform(size=n_apa + 1) + 1
                ws[n_apa] += -1
                ws = ws / np.sum(ws)
                theta_win_mat = self.init_theta_win(n_apa)
                if self.verbose:
                    logger.debug(
                        f'''
                        init alpha: {theta_win_mat.res_alpha_arr.values.tolist()},
                        init beta: {theta_win_mat.res_beta_arr.values.tolist()}
                        ''')

                res_list[i] = \
                    EM.em_algo(ws,
                               theta_win_mat,
                               self.n_frag,
                               self.pre_alpha_arr,
                               self.pre_beta_arr,
                               self.all_theta,
                               self.init_beta_arr,
                               self.all_win_log_lik_mat_list,
                               self.pre_init_alpha_arr,
                               self.min_pa_gap,
                               self.unif_log_lik,
                               self.min_theta,
                               self.L,
                               verbose=verbose)
            except Exception as e:
                traceback.print_exc()
                logger.debug(f'Error found in {i+1} trial. Next')
                res_list[i] = \
                    dotdict({
                        'ws': None,
                        'alpha_arr': None,
                        'beta_arr': None,
                        'lb_arr': None,
                        'bic': None
                    })
            if not isinstance(res_list[i].ws, np.ndarray):
                continue
            # 20200728 update
            # print(np.diff(res_list[i].alpha_arr))
            if any(np.diff(res_list[i].alpha_arr) < self.min_pa_gap):
                continue

            lb_arr[i] = res_list[i].lb_arr[-1]
            bic_arr[i] = res_list[i].bic
            if verbose:
                logger.debug(
                    f'''
                        K = {n_apa}, {i+1}_trial, n_trial_{n_trial}:
                        ws: {np.round(res_list[i].ws, decimals=3)}
                        alpha: {res_list[i].alpha_arr.values.tolist()}
                        beta: {res_list[i].beta_arr.values.tolist()}
                        lb = {res_list[i].lb_arr[-1]}
                        bic = {res_list[i].bic}
                        '''
                )


        min_ind = np.argmin(bic_arr)
        res = res_list[min_ind]

        if isinstance(res.ws, np.ndarray):
            if verbose:
                logger.debug(
                    f'''
                        K = {n_apa}, final result:
                        ws: {np.round(res.ws, decimals=3)}
                        alpha: {res.alpha_arr.values.tolist()}
                        beta: {res.beta_arr.values.tolist()}
                        lb = {res.lb_arr[-1]}
                        bic = {res.bic}
                        '''
                )
        else:
            logger.debug(
                f'''
                    No results available for K = {n_apa}
                '''
            )

        return res

    def inference(self):
        u'''
        return the given init_theta_win_mat for model, but no need
        '''
        if self.fixed_inference_flag and not self.pre_init_theta_win_mat:

            alpha_arr = replace_with_closest(
                self.all_theta,
                self.pre_init_theta_win_mat[:, 0]
            )
            beta_arr = replace_with_closest(
                self.all_theta,
                self.pre_init_theta_win_mat[:, 1]
            )
            res = EM.fixed_inference(alpha_arr,
                                     self.n_frag,
                                     self.nround,
                                     self.all_theta,
                                     self.all_win_log_lik_mat_list,
                                     self.init_beta_arr,
                                     self.unif_log_lik,
                                     verbose=self.verbose
                                     )
            res.alpha_arr = self.pre_init_theta_win_mat[:, 0]
            res.beta_arr = self.pre_init_theta_win_mat[:, 0]
            return res

        res_list = []
        lb_arr = np.repeat(-np.Inf, self.n_max_apa)
        bic_arr = np.repeat(np.Inf, self.n_max_apa)

        for i in range(self.n_max_apa, self.n_min_apa - 1, -1):
            res_tmp = self.em_optim(n_apa=i, verbose=self.verbose)
            res_list.append(res_tmp)
            if not isinstance(res_tmp.ws, np.ndarray):
                continue
            lb_arr[i - 1] = res_tmp.lb_arr[-1]
            bic_arr[i - 1] = res_tmp.bic

        res_list = res_list[::-1]
        min_ind = np.argmin(bic_arr)
        res = res_list[min_ind]


        if not res:
            logger.warning(
                f'''
                Inference failed. No results available.
                '''
                )

        res_2 = self.rm_component(res,
                                  self.min_ws,
                                  self.pre_alpha_arr,
                                  self.pre_init_alpha_arr,
                                  self.n_frag,
                                  self.all_theta,
                                  self.all_win_log_lik_mat_list,
                                  self.init_beta_arr,
                                  self.unif_log_lik,
                                  verbose=self.verbose)
        if self.verbose:
            logger.debug(
                f'''
                    Final result:
                    ws: {np.round(res_2.ws, decimals=3)}
                    alpha: {res_2.alpha_arr.values.tolist()}
                    beta: {res_2.beta_arr.values.tolist()}
                    lb = {res_2.lb_arr[-1]}
                    bic = {res_2.bic}
                    label= {res_2.label}
                '''
            )

        return res_2

    @staticmethod
    def rm_component(res,
                     min_ws,
                     pre_alpha_arr,
                     pre_init_alpha_arr,
                     n_frag,
                     all_theta,
                     all_win_log_lik_mat_list,
                     init_beta_arr,
                     unif_log_lik,
                     nround=200,
                     verbose=False):
        K = len(res.alpha_arr)
        if pre_alpha_arr:
            pre_flag = np.isin(
                res.alpha_arr, pre_init_alpha_arr
            )
            kept_inds = np.invert(np.logical_and(
                res.ws[:K] < min_ws,
                np.invert(pre_flag)
            ))
        else:
            kept_inds = res.ws[:K] >= min_ws

        if np.sum(kept_inds) == K or np.sum(kept_inds) == 0:
            log_zmat = np.zeros(
                (n_frag, K + 1)
                )
            for k in range(K + 1):
                log_zmat = EM.cal_z_k(
                    all_win_log_lik_mat_list,
                    res.alpha_arr,
                    res.beta_arr,
                    res.ws,
                    k,
                    log_zmat,
                    all_theta,
                    init_beta_arr,
                    unif_log_lik
                )
            Z = EM.norm_z(log_zmat)
            label = log_zmat.argmax(axis=1)
            res.label = label
            return res

        else:
            alpha_arr = res.alpha_arr[kept_inds].reset_index(drop=True)
            beta_arr = res.beta_arr[kept_inds].reset_index(drop=True)
            res = EM.fixed_inference(alpha_arr,
                                     beta_arr,
                                     n_frag,
                                     nround,
                                     all_theta,
                                     all_win_log_lik_mat_list,
                                     init_beta_arr,
                                     unif_log_lik,
                                     verbose=verbose
                                     )
            return res
