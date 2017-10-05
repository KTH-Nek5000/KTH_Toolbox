!> @file rprm_block.f
!! @ingroup runparam
!! @brief Block data to initialise common block for runtime parameter module
!! @details Following Nek5000 standard I keep block data in seaprate file.
!! @author Adam Peplinski
!! @date Sep 28, 2017
!=======================================================================
      block data rprm_common_init
      include 'RPRMD'

      data rprm_ifinit /.false./
      data rprm_pid0 /0/
      data rprm_par_num /0/
      data rprm_par_mpos /0/
      data rprm_par_id /rprm_id_size*-1/
      data rprm_par_name /rprm_id_max*rprm_blname/
      data rprm_parv_int /rprm_id_max*0/
      data rprm_parv_real /rprm_id_max*0.0/
      data rprm_parv_log /rprm_id_max*.false./
      data rprm_parv_str /rprm_id_max*rprm_blname/

      end
