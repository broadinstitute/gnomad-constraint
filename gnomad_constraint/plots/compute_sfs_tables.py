{
  "cells": [
    {
      "cell_type": "markdown",
      "metadata": {},
      "source": [
        "# Site Frequency Spectrum (SFS) Computation\n",
        "\n",
        "This notebook computes Site Frequency Spectrum tables for each combination of mutation type and functional effect across different downsampling levels.\n",
        "\n",
        "## Output Format\n",
        "For each combination of functional category (syn, mis, lof) and mutation rate category (transversion, transition, CpG), we generate a table named:\n",
        "`SFS_{functional_category}_{mutation_rate_category}.txt.gz`\n",
        "\n",
        "Each table contains:\n",
        "- Rows: Different downsampling levels (AC_ds10, AC_ds100, AC_ds500, etc.)\n",
        "- Columns: Allele count bins (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10-100, 100-1000, 1000-10000, >10000)\n",
        "- Values: Number of variants with that allele count in that downsampling\n"
      ]
    },
    {
      "cell_type": "code",
      "execution_count": 17,
      "metadata": {},
      "outputs": [
        {
          "name": "stderr",
          "output_type": "stream",
          "text": [
            "WARNING (backend.service_backend 192): Hail has already been initialized. If this call was intended to change configuration, close the session with hl.stop() first.\n",
            "/Users/jgoodric/miniconda3/envs/loftee/lib/python3.11/site-packages/hail/backend/service_backend.py:563: UserWarning:\n",
            "\n",
            "Modifying the requester pays project or buckets at runtime using flags is deprecated. Expect this behavior to become unsupported soon.\n",
            "\n"
          ]
        },
        {
          "ename": "AssertionError",
          "evalue": "",
          "output_type": "error",
          "traceback": [
            "\u001b[31m---------------------------------------------------------------------------\u001b[39m",
            "\u001b[31mAssertionError\u001b[39m                            Traceback (most recent call last)",
            "\u001b[36mCell\u001b[39m\u001b[36m \u001b[39m\u001b[32mIn[17]\u001b[39m\u001b[32m, line 9\u001b[39m\n\u001b[32m      6\u001b[39m \u001b[38;5;28;01mimport\u001b[39;00m\u001b[38;5;250m \u001b[39m\u001b[34;01mos\u001b[39;00m\n\u001b[32m      8\u001b[39m \u001b[38;5;66;03m# Initialize Hail\u001b[39;00m\n\u001b[32m----> \u001b[39m\u001b[32m9\u001b[39m \u001b[43mhl\u001b[49m\u001b[43m.\u001b[49m\u001b[43minit\u001b[49m\u001b[43m(\u001b[49m\u001b[43mlog\u001b[49m\u001b[43m=\u001b[49m\u001b[33;43m'\u001b[39;49m\u001b[33;43m/tmp/sfs_computation.log\u001b[39;49m\u001b[33;43m'\u001b[39;49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mtmp_dir\u001b[49m\u001b[43m=\u001b[49m\u001b[33;43m'\u001b[39;49m\u001b[33;43mgs://gnomad-tmp-4day\u001b[39;49m\u001b[33;43m'\u001b[39;49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mbackend\u001b[49m\u001b[43m=\u001b[49m\u001b[33;43m\"\u001b[39;49m\u001b[33;43mbatch\u001b[39;49m\u001b[33;43m\"\u001b[39;49m\u001b[43m)\u001b[49m\n\u001b[32m     11\u001b[39m \u001b[38;5;28mprint\u001b[39m(\u001b[33m\"\u001b[39m\u001b[33mHail initialized successfully\u001b[39m\u001b[33m\"\u001b[39m)\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/decorator.py:235\u001b[39m, in \u001b[36mdecorate.<locals>.fun\u001b[39m\u001b[34m(*args, **kw)\u001b[39m\n\u001b[32m    233\u001b[39m \u001b[38;5;28;01mif\u001b[39;00m \u001b[38;5;129;01mnot\u001b[39;00m kwsyntax:\n\u001b[32m    234\u001b[39m     args, kw = fix(args, kw, sig)\n\u001b[32m--> \u001b[39m\u001b[32m235\u001b[39m \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[43mcaller\u001b[49m\u001b[43m(\u001b[49m\u001b[43mfunc\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43m*\u001b[49m\u001b[43m(\u001b[49m\u001b[43mextras\u001b[49m\u001b[43m \u001b[49m\u001b[43m+\u001b[49m\u001b[43m \u001b[49m\u001b[43margs\u001b[49m\u001b[43m)\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43m*\u001b[49m\u001b[43m*\u001b[49m\u001b[43mkw\u001b[49m\u001b[43m)\u001b[49m\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/hail/typecheck/check.py:585\u001b[39m, in \u001b[36m_make_dec.<locals>.wrapper\u001b[39m\u001b[34m(__original_func, *args, **kwargs)\u001b[39m\n\u001b[32m    582\u001b[39m \u001b[38;5;129m@decorator\u001b[39m\n\u001b[32m    583\u001b[39m \u001b[38;5;28;01mdef\u001b[39;00m\u001b[38;5;250m \u001b[39m\u001b[34mwrapper\u001b[39m(__original_func: Callable[..., T], *args, **kwargs) -> T:\n\u001b[32m    584\u001b[39m     args_, kwargs_ = check_all(__original_func, args, kwargs, checkers, is_method=is_method)\n\u001b[32m--> \u001b[39m\u001b[32m585\u001b[39m     \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[43m__original_func\u001b[49m\u001b[43m(\u001b[49m\u001b[43m*\u001b[49m\u001b[43margs_\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43m*\u001b[49m\u001b[43m*\u001b[49m\u001b[43mkwargs_\u001b[49m\u001b[43m)\u001b[49m\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/hail/context.py:371\u001b[39m, in \u001b[36minit\u001b[39m\u001b[34m(sc, app_name, master, local, log, quiet, append, min_block_size, branching_factor, tmp_dir, default_reference, idempotent, global_seed, spark_conf, skip_logging_configuration, local_tmpdir, _optimizer_iterations, backend, driver_cores, driver_memory, worker_cores, worker_memory, gcs_requester_pays_configuration, regions, gcs_bucket_allow_list, copy_spark_log_on_error)\u001b[39m\n\u001b[32m    368\u001b[39m     backend = \u001b[33m'\u001b[39m\u001b[33mbatch\u001b[39m\u001b[33m'\u001b[39m\n\u001b[32m    370\u001b[39m \u001b[38;5;28;01mif\u001b[39;00m backend == \u001b[33m'\u001b[39m\u001b[33mbatch\u001b[39m\u001b[33m'\u001b[39m:\n\u001b[32m--> \u001b[39m\u001b[32m371\u001b[39m     \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[43mhail_event_loop\u001b[49m\u001b[43m(\u001b[49m\u001b[43m)\u001b[49m\u001b[43m.\u001b[49m\u001b[43mrun_until_complete\u001b[49m\u001b[43m(\u001b[49m\n\u001b[32m    372\u001b[39m \u001b[43m        \u001b[49m\u001b[43minit_batch\u001b[49m\u001b[43m(\u001b[49m\n\u001b[32m    373\u001b[39m \u001b[43m            \u001b[49m\u001b[43mlog\u001b[49m\u001b[43m=\u001b[49m\u001b[43mlog\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m    374\u001b[39m \u001b[43m            \u001b[49m\u001b[43mquiet\u001b[49m\u001b[43m=\u001b[49m\u001b[43mquiet\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m    375\u001b[39m \u001b[43m            \u001b[49m\u001b[43mappend\u001b[49m\u001b[43m=\u001b[49m\u001b[43mappend\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m    376\u001b[39m \u001b[43m            \u001b[49m\u001b[43mtmpdir\u001b[49m\u001b[43m=\u001b[49m\u001b[43mtmp_dir\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m    377\u001b[39m \u001b[43m            \u001b[49m\u001b[43mlocal_tmpdir\u001b[49m\u001b[43m=\u001b[49m\u001b[43mlocal_tmpdir\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m    378\u001b[39m \u001b[43m            \u001b[49m\u001b[43mdefault_reference\u001b[49m\u001b[43m=\u001b[49m\u001b[43mdefault_reference\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m    379\u001b[39m \u001b[43m            \u001b[49m\u001b[43mglobal_seed\u001b[49m\u001b[43m=\u001b[49m\u001b[43mglobal_seed\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m    380\u001b[39m \u001b[43m            \u001b[49m\u001b[43mdriver_cores\u001b[49m\u001b[43m=\u001b[49m\u001b[43mdriver_cores\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m    381\u001b[39m \u001b[43m            \u001b[49m\u001b[43mdriver_memory\u001b[49m\u001b[43m=\u001b[49m\u001b[43mdriver_memory\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m    382\u001b[39m \u001b[43m            \u001b[49m\u001b[43mworker_cores\u001b[49m\u001b[43m=\u001b[49m\u001b[43mworker_cores\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m    383\u001b[39m \u001b[43m            \u001b[49m\u001b[43mworker_memory\u001b[49m\u001b[43m=\u001b[49m\u001b[43mworker_memory\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m    384\u001b[39m \u001b[43m            \u001b[49m\u001b[43mname_prefix\u001b[49m\u001b[43m=\u001b[49m\u001b[43mapp_name\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m    385\u001b[39m \u001b[43m            \u001b[49m\u001b[43mgcs_requester_pays_configuration\u001b[49m\u001b[43m=\u001b[49m\u001b[43mgcs_requester_pays_configuration\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m    386\u001b[39m \u001b[43m            \u001b[49m\u001b[43mregions\u001b[49m\u001b[43m=\u001b[49m\u001b[43mregions\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m    387\u001b[39m \u001b[43m            \u001b[49m\u001b[43mgcs_bucket_allow_list\u001b[49m\u001b[43m=\u001b[49m\u001b[43mgcs_bucket_allow_list\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m    388\u001b[39m \u001b[43m        \u001b[49m\u001b[43m)\u001b[49m\n\u001b[32m    389\u001b[39m \u001b[43m    \u001b[49m\u001b[43m)\u001b[49m\n\u001b[32m    390\u001b[39m \u001b[38;5;28;01mif\u001b[39;00m backend == \u001b[33m'\u001b[39m\u001b[33mspark\u001b[39m\u001b[33m'\u001b[39m:\n\u001b[32m    391\u001b[39m     \u001b[38;5;28;01mreturn\u001b[39;00m init_spark(\n\u001b[32m    392\u001b[39m         sc=sc,\n\u001b[32m    393\u001b[39m         app_name=app_name,\n\u001b[32m   (...)\u001b[39m\u001b[32m    410\u001b[39m         copy_log_on_error=copy_spark_log_on_error,\n\u001b[32m    411\u001b[39m     )\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/nest_asyncio.py:98\u001b[39m, in \u001b[36m_patch_loop.<locals>.run_until_complete\u001b[39m\u001b[34m(self, future)\u001b[39m\n\u001b[32m     95\u001b[39m \u001b[38;5;28;01mif\u001b[39;00m \u001b[38;5;129;01mnot\u001b[39;00m f.done():\n\u001b[32m     96\u001b[39m     \u001b[38;5;28;01mraise\u001b[39;00m \u001b[38;5;167;01mRuntimeError\u001b[39;00m(\n\u001b[32m     97\u001b[39m         \u001b[33m'\u001b[39m\u001b[33mEvent loop stopped before Future completed.\u001b[39m\u001b[33m'\u001b[39m)\n\u001b[32m---> \u001b[39m\u001b[32m98\u001b[39m \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[43mf\u001b[49m\u001b[43m.\u001b[49m\u001b[43mresult\u001b[49m\u001b[43m(\u001b[49m\u001b[43m)\u001b[49m\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/asyncio/futures.py:203\u001b[39m, in \u001b[36mFuture.result\u001b[39m\u001b[34m(self)\u001b[39m\n\u001b[32m    201\u001b[39m \u001b[38;5;28mself\u001b[39m.__log_traceback = \u001b[38;5;28;01mFalse\u001b[39;00m\n\u001b[32m    202\u001b[39m \u001b[38;5;28;01mif\u001b[39;00m \u001b[38;5;28mself\u001b[39m._exception \u001b[38;5;129;01mis\u001b[39;00m \u001b[38;5;129;01mnot\u001b[39;00m \u001b[38;5;28;01mNone\u001b[39;00m:\n\u001b[32m--> \u001b[39m\u001b[32m203\u001b[39m     \u001b[38;5;28;01mraise\u001b[39;00m \u001b[38;5;28mself\u001b[39m._exception.with_traceback(\u001b[38;5;28mself\u001b[39m._exception_tb)\n\u001b[32m    204\u001b[39m \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[38;5;28mself\u001b[39m._result\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/asyncio/tasks.py:267\u001b[39m, in \u001b[36mTask.__step\u001b[39m\u001b[34m(***failed resolving arguments***)\u001b[39m\n\u001b[32m    263\u001b[39m \u001b[38;5;28;01mtry\u001b[39;00m:\n\u001b[32m    264\u001b[39m     \u001b[38;5;28;01mif\u001b[39;00m exc \u001b[38;5;129;01mis\u001b[39;00m \u001b[38;5;28;01mNone\u001b[39;00m:\n\u001b[32m    265\u001b[39m         \u001b[38;5;66;03m# We use the `send` method directly, because coroutines\u001b[39;00m\n\u001b[32m    266\u001b[39m         \u001b[38;5;66;03m# don't have `__iter__` and `__next__` methods.\u001b[39;00m\n\u001b[32m--> \u001b[39m\u001b[32m267\u001b[39m         result = coro.send(\u001b[38;5;28;01mNone\u001b[39;00m)\n\u001b[32m    268\u001b[39m     \u001b[38;5;28;01melse\u001b[39;00m:\n\u001b[32m    269\u001b[39m         result = coro.throw(exc)\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/hail/context.py:579\u001b[39m, in \u001b[36minit_batch\u001b[39m\u001b[34m(billing_project, remote_tmpdir, log, quiet, append, tmpdir, local_tmpdir, default_reference, global_seed, disable_progress_bar, driver_cores, driver_memory, worker_cores, worker_memory, name_prefix, token, gcs_requester_pays_configuration, regions, gcs_bucket_allow_list)\u001b[39m\n\u001b[32m    576\u001b[39m     tmpdir = backend.remote_tmpdir + \u001b[33m'\u001b[39m\u001b[33mtmp/hail/\u001b[39m\u001b[33m'\u001b[39m + secret_alnum_string()\n\u001b[32m    577\u001b[39m local_tmpdir = _get_local_tmpdir(local_tmpdir)\n\u001b[32m--> \u001b[39m\u001b[32m579\u001b[39m \u001b[43mHailContext\u001b[49m\u001b[43m.\u001b[49m\u001b[43mcreate\u001b[49m\u001b[43m(\u001b[49m\u001b[43mlog\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mquiet\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mappend\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mtmpdir\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mlocal_tmpdir\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mdefault_reference\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mglobal_seed\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mbackend\u001b[49m\u001b[43m)\u001b[49m\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/hail/context.py:76\u001b[39m, in \u001b[36mHailContext.create\u001b[39m\u001b[34m(log, quiet, append, tmpdir, local_tmpdir, default_reference, global_seed, backend)\u001b[39m\n\u001b[32m     65\u001b[39m \u001b[38;5;129m@staticmethod\u001b[39m\n\u001b[32m     66\u001b[39m \u001b[38;5;28;01mdef\u001b[39;00m\u001b[38;5;250m \u001b[39m\u001b[34mcreate\u001b[39m(\n\u001b[32m     67\u001b[39m     log: \u001b[38;5;28mstr\u001b[39m,\n\u001b[32m   (...)\u001b[39m\u001b[32m     74\u001b[39m     backend: Backend,\n\u001b[32m     75\u001b[39m ):\n\u001b[32m---> \u001b[39m\u001b[32m76\u001b[39m     hc = \u001b[43mHailContext\u001b[49m\u001b[43m(\u001b[49m\n\u001b[32m     77\u001b[39m \u001b[43m        \u001b[49m\u001b[43mlog\u001b[49m\u001b[43m=\u001b[49m\u001b[43mlog\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m     78\u001b[39m \u001b[43m        \u001b[49m\u001b[43mquiet\u001b[49m\u001b[43m=\u001b[49m\u001b[43mquiet\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m     79\u001b[39m \u001b[43m        \u001b[49m\u001b[43mappend\u001b[49m\u001b[43m=\u001b[49m\u001b[43mappend\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m     80\u001b[39m \u001b[43m        \u001b[49m\u001b[43mtmpdir\u001b[49m\u001b[43m=\u001b[49m\u001b[43mtmpdir\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m     81\u001b[39m \u001b[43m        \u001b[49m\u001b[43mlocal_tmpdir\u001b[49m\u001b[43m=\u001b[49m\u001b[43mlocal_tmpdir\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m     82\u001b[39m \u001b[43m        \u001b[49m\u001b[43mglobal_seed\u001b[49m\u001b[43m=\u001b[49m\u001b[43mglobal_seed\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m     83\u001b[39m \u001b[43m        \u001b[49m\u001b[43mbackend\u001b[49m\u001b[43m=\u001b[49m\u001b[43mbackend\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m     84\u001b[39m \u001b[43m    \u001b[49m\u001b[43m)\u001b[49m\n\u001b[32m     85\u001b[39m     hc.initialize_references(default_reference)\n\u001b[32m     86\u001b[39m     \u001b[38;5;28;01mreturn\u001b[39;00m hc\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/decorator.py:235\u001b[39m, in \u001b[36mdecorate.<locals>.fun\u001b[39m\u001b[34m(*args, **kw)\u001b[39m\n\u001b[32m    233\u001b[39m \u001b[38;5;28;01mif\u001b[39;00m \u001b[38;5;129;01mnot\u001b[39;00m kwsyntax:\n\u001b[32m    234\u001b[39m     args, kw = fix(args, kw, sig)\n\u001b[32m--> \u001b[39m\u001b[32m235\u001b[39m \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[43mcaller\u001b[49m\u001b[43m(\u001b[49m\u001b[43mfunc\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43m*\u001b[49m\u001b[43m(\u001b[49m\u001b[43mextras\u001b[49m\u001b[43m \u001b[49m\u001b[43m+\u001b[49m\u001b[43m \u001b[49m\u001b[43margs\u001b[49m\u001b[43m)\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43m*\u001b[49m\u001b[43m*\u001b[49m\u001b[43mkw\u001b[49m\u001b[43m)\u001b[49m\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/hail/typecheck/check.py:585\u001b[39m, in \u001b[36m_make_dec.<locals>.wrapper\u001b[39m\u001b[34m(__original_func, *args, **kwargs)\u001b[39m\n\u001b[32m    582\u001b[39m \u001b[38;5;129m@decorator\u001b[39m\n\u001b[32m    583\u001b[39m \u001b[38;5;28;01mdef\u001b[39;00m\u001b[38;5;250m \u001b[39m\u001b[34mwrapper\u001b[39m(__original_func: Callable[..., T], *args, **kwargs) -> T:\n\u001b[32m    584\u001b[39m     args_, kwargs_ = check_all(__original_func, args, kwargs, checkers, is_method=is_method)\n\u001b[32m--> \u001b[39m\u001b[32m585\u001b[39m     \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[43m__original_func\u001b[49m\u001b[43m(\u001b[49m\u001b[43m*\u001b[49m\u001b[43margs_\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43m*\u001b[49m\u001b[43m*\u001b[49m\u001b[43mkwargs_\u001b[49m\u001b[43m)\u001b[49m\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/hail/context.py:92\u001b[39m, in \u001b[36mHailContext.__init__\u001b[39m\u001b[34m(self, log, quiet, append, tmpdir, local_tmpdir, global_seed, backend)\u001b[39m\n\u001b[32m     88\u001b[39m \u001b[38;5;129m@typecheck_method\u001b[39m(\n\u001b[32m     89\u001b[39m     log=\u001b[38;5;28mstr\u001b[39m, quiet=\u001b[38;5;28mbool\u001b[39m, append=\u001b[38;5;28mbool\u001b[39m, tmpdir=\u001b[38;5;28mstr\u001b[39m, local_tmpdir=\u001b[38;5;28mstr\u001b[39m, global_seed=nullable(\u001b[38;5;28mint\u001b[39m), backend=Backend\n\u001b[32m     90\u001b[39m )\n\u001b[32m     91\u001b[39m \u001b[38;5;28;01mdef\u001b[39;00m\u001b[38;5;250m \u001b[39m\u001b[34m__init__\u001b[39m(\u001b[38;5;28mself\u001b[39m, log, quiet, append, tmpdir, local_tmpdir, global_seed, backend):\n\u001b[32m---> \u001b[39m\u001b[32m92\u001b[39m     \u001b[38;5;28;01massert\u001b[39;00m \u001b[38;5;129;01mnot\u001b[39;00m Env._hc\n\u001b[32m     94\u001b[39m     \u001b[38;5;28mself\u001b[39m._log = log\n\u001b[32m     96\u001b[39m     \u001b[38;5;28mself\u001b[39m._tmpdir = tmpdir\n",
            "\u001b[31mAssertionError\u001b[39m: "
          ]
        }
      ],
      "source": [
        "import hail as hl\n",
        "import pandas as pd\n",
        "import numpy as np\n",
        "from typing import Dict, List, Tuple\n",
        "import gzip\n",
        "import os\n",
        "\n",
        "# Initialize Hail\n",
        "hl.init(log='/tmp/sfs_computation.log', tmp_dir='gs://gnomad-tmp-4day', backend=\"batch\")\n",
        "\n",
        "print(\"Hail initialized successfully\")\n"
      ]
    },
    {
      "cell_type": "markdown",
      "metadata": {},
      "source": [
        "## Load and Prepare Data\n",
        "\n",
        "Load the constraint data and prepare it for SFS computation.\n"
      ]
    },
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {},
      "outputs": [
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\"></pre>\n"
            ],
            "text/plain": [
              "\u001b[?25l"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "application/vnd.jupyter.widget-view+json": {
              "model_id": "1dc85a7644fb40d688560ac1a686cb75",
              "version_major": 2,
              "version_minor": 0
            },
            "text/plain": [
              "Output()"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\">\n",
              "</pre>\n"
            ],
            "text/plain": [
              "\n",
              "\u001b[?25h"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\"></pre>\n"
            ],
            "text/plain": [
              "\u001b[?25l"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "application/vnd.jupyter.widget-view+json": {
              "model_id": "214a694b396e406984f76e51d94bb8f7",
              "version_major": 2,
              "version_minor": 0
            },
            "text/plain": [
              "Output()"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\">\n",
              "</pre>\n"
            ],
            "text/plain": [
              "\n",
              "\u001b[?25h"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\"></pre>\n"
            ],
            "text/plain": [
              "\u001b[?25l"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "application/vnd.jupyter.widget-view+json": {
              "model_id": "bbfed77d25ee4c68b00dbb0c6b433405",
              "version_major": 2,
              "version_minor": 0
            },
            "text/plain": [
              "Output()"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\">\n",
              "</pre>\n"
            ],
            "text/plain": [
              "\n",
              "\u001b[?25h"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\"></pre>\n"
            ],
            "text/plain": [
              "\u001b[?25l"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "application/vnd.jupyter.widget-view+json": {
              "model_id": "a9a17be2333147b9871f48fbd9cdc2cc",
              "version_major": 2,
              "version_minor": 0
            },
            "text/plain": [
              "Output()"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\">\n",
              "</pre>\n"
            ],
            "text/plain": [
              "\n",
              "\u001b[?25h"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "name": "stdout",
          "output_type": "stream",
          "text": [
            "\n",
            "Table schema:\n",
            "----------------------------------------\n",
            "Global fields:\n",
            "    'calculate_mu_globals': struct {\n",
            "        freq_meta: array<dict<str, str>>, \n",
            "        ac_cutoff: int32, \n",
            "        min_cov: int32, \n",
            "        max_cov: int32, \n",
            "        gerp_lower_cutoff: float64, \n",
            "        gerp_upper_cutoff: float64, \n",
            "        genetic_ancestry_groups: array<str>, \n",
            "        downsampling_level: int32, \n",
            "        downsampling_idx: int32, \n",
            "        most_severe_consequence: array<str>\n",
            "    } \n",
            "    'build_models_globals': struct {\n",
            "        synonymous_transcript_filter_field: str, \n",
            "        low_cov_cutoff: int32, \n",
            "        high_cov_cutoff: int32, \n",
            "        upper_cov_cutoff: int32, \n",
            "        skip_coverage_model: bool\n",
            "    } \n",
            "    'apply_models_globals': struct {\n",
            "        low_cov_cutoff: int32, \n",
            "        high_cov_cutoff: int32, \n",
            "        skip_coverage_model: bool\n",
            "    } \n",
            "    'exomes_freq_meta': array<dict<str, str>> \n",
            "    'genetic_ancestry_groups': array<str> \n",
            "    'downsamplings': array<int32> \n",
            "    'max_af': float64 \n",
            "----------------------------------------\n",
            "Row fields:\n",
            "    'locus': locus<GRCh38> \n",
            "    'alleles': array<str> \n",
            "    'mutation_type': str \n",
            "    'most_severe_consequence': struct {\n",
            "        transcript_id: str, \n",
            "        gene_id: str, \n",
            "        gene_symbol: str, \n",
            "        biotype: str, \n",
            "        most_severe_consequence: str, \n",
            "        mane_select: str, \n",
            "        canonical: int32, \n",
            "        lof: str\n",
            "    } \n",
            "    'canonical_most_severe_consequence': struct {\n",
            "        transcript_id: str, \n",
            "        gene_id: str, \n",
            "        gene_symbol: str, \n",
            "        biotype: str, \n",
            "        most_severe_consequence: str, \n",
            "        mane_select: str, \n",
            "        canonical: int32, \n",
            "        lof: str\n",
            "    } \n",
            "    'mane_select_most_severe_consequence': struct {\n",
            "        transcript_id: str, \n",
            "        gene_id: str, \n",
            "        gene_symbol: str, \n",
            "        biotype: str, \n",
            "        most_severe_consequence: str, \n",
            "        mane_select: str, \n",
            "        canonical: int32, \n",
            "        lof: str\n",
            "    } \n",
            "    'AC': int32 \n",
            "    'AF': float64 \n",
            "    'ac_group': str \n",
            "    'downsampling': str \n",
            "----------------------------------------\n",
            "Key: ['locus', 'alleles']\n",
            "----------------------------------------\n"
          ]
        },
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\"></pre>\n"
            ],
            "text/plain": [
              "\u001b[?25l"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "application/vnd.jupyter.widget-view+json": {
              "model_id": "792f056043bf4b0fb4c47054fd753ec8",
              "version_major": 2,
              "version_minor": 0
            },
            "text/plain": [
              "Output()"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\">\n",
              "</pre>\n"
            ],
            "text/plain": [
              "\n",
              "\u001b[?25h"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "text/html": [
              "<table><thead><tr><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;\"></div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;\"></div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;\"></div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"8\"><div style=\"text-align: left;\"></div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"8\"><div style=\"text-align: left;\"></div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"8\"><div style=\"text-align: left;\"></div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;\"></div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;\"></div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;\"></div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;\"></div></td></tr><tr><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;\"></div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;\"></div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;\"></div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"8\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">most_severe_consequence</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"8\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">canonical_most_severe_consequence</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"8\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">mane_select_most_severe_consequence</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;\"></div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;\"></div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;\"></div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;\"></div></td></tr><tr><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">locus</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">alleles</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">mutation_type</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">transcript_id</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">gene_id</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">gene_symbol</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">biotype</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">most_severe_consequence</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">mane_select</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">canonical</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">lof</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">transcript_id</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">gene_id</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">gene_symbol</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">biotype</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">most_severe_consequence</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">mane_select</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">canonical</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">lof</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">transcript_id</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">gene_id</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">gene_symbol</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">biotype</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">most_severe_consequence</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">mane_select</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">canonical</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">lof</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">AC</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">AF</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">ac_group</div></td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \" colspan=\"1\"><div style=\"text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px\">downsampling</div></td></tr><tr><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">locus&lt;GRCh38&gt;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">array&lt;str&gt;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">int32</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">int32</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">int32</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">int32</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">float64</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;\">str</td></tr>\n",
              "</thead><tbody><tr><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">chr1:69402</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">[&quot;G&quot;,&quot;T&quot;]</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;transversion&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0.00e+00</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;0&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;AC_ds10&quot;</td></tr>\n",
              "<tr><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">chr1:69402</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">[&quot;G&quot;,&quot;T&quot;]</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;transversion&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0.00e+00</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;0&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;AC_ds100&quot;</td></tr>\n",
              "<tr><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">chr1:69402</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">[&quot;G&quot;,&quot;T&quot;]</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;transversion&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0.00e+00</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;0&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;AC_ds500&quot;</td></tr>\n",
              "<tr><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">chr1:69402</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">[&quot;G&quot;,&quot;T&quot;]</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;transversion&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0.00e+00</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;0&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;AC_ds1000&quot;</td></tr>\n",
              "<tr><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">chr1:69402</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">[&quot;G&quot;,&quot;T&quot;]</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;transversion&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0.00e+00</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;0&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;AC_ds2000&quot;</td></tr>\n",
              "<tr><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">chr1:69402</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">[&quot;G&quot;,&quot;T&quot;]</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;transversion&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0.00e+00</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;0&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;AC_ds2884&quot;</td></tr>\n",
              "<tr><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">chr1:69402</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">[&quot;G&quot;,&quot;T&quot;]</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;transversion&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0.00e+00</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;0&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;AC_ds5000&quot;</td></tr>\n",
              "<tr><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">chr1:69402</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">[&quot;G&quot;,&quot;T&quot;]</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;transversion&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0.00e+00</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;0&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;AC_ds10000&quot;</td></tr>\n",
              "<tr><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">chr1:69402</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">[&quot;G&quot;,&quot;T&quot;]</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;transversion&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0.00e+00</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;0&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;AC_ds13068&quot;</td></tr>\n",
              "<tr><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">chr1:69402</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">[&quot;G&quot;,&quot;T&quot;]</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;transversion&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENST00000641515&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;ENSG00000186092&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;OR4F5&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;protein_coding&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;missense_variant&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;NM_001005484.2&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">1</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">NA</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">0.00e+00</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;0&quot;</td><td style=\"white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; \">&quot;AC_ds16740&quot;</td></tr>\n",
              "</tbody></table><p style=\"background: #fdd; padding: 0.4em;\">showing top 10 rows</p>\n"
            ],
            "text/plain": [
              "+---------------+------------+----------------+\n",
              "| locus         | alleles    | mutation_type  |\n",
              "+---------------+------------+----------------+\n",
              "| locus<GRCh38> | array<str> | str            |\n",
              "+---------------+------------+----------------+\n",
              "| chr1:69402    | [\"G\",\"T\"]  | \"transversion\" |\n",
              "| chr1:69402    | [\"G\",\"T\"]  | \"transversion\" |\n",
              "| chr1:69402    | [\"G\",\"T\"]  | \"transversion\" |\n",
              "| chr1:69402    | [\"G\",\"T\"]  | \"transversion\" |\n",
              "| chr1:69402    | [\"G\",\"T\"]  | \"transversion\" |\n",
              "| chr1:69402    | [\"G\",\"T\"]  | \"transversion\" |\n",
              "| chr1:69402    | [\"G\",\"T\"]  | \"transversion\" |\n",
              "| chr1:69402    | [\"G\",\"T\"]  | \"transversion\" |\n",
              "| chr1:69402    | [\"G\",\"T\"]  | \"transversion\" |\n",
              "| chr1:69402    | [\"G\",\"T\"]  | \"transversion\" |\n",
              "+---------------+------------+----------------+\n",
              "\n",
              "+---------------------------------------+---------------------------------+\n",
              "| most_severe_consequence.transcript_id | most_severe_consequence.gene_id |\n",
              "+---------------------------------------+---------------------------------+\n",
              "| str                                   | str                             |\n",
              "+---------------------------------------+---------------------------------+\n",
              "| \"ENST00000641515\"                     | \"ENSG00000186092\"               |\n",
              "| \"ENST00000641515\"                     | \"ENSG00000186092\"               |\n",
              "| \"ENST00000641515\"                     | \"ENSG00000186092\"               |\n",
              "| \"ENST00000641515\"                     | \"ENSG00000186092\"               |\n",
              "| \"ENST00000641515\"                     | \"ENSG00000186092\"               |\n",
              "| \"ENST00000641515\"                     | \"ENSG00000186092\"               |\n",
              "| \"ENST00000641515\"                     | \"ENSG00000186092\"               |\n",
              "| \"ENST00000641515\"                     | \"ENSG00000186092\"               |\n",
              "| \"ENST00000641515\"                     | \"ENSG00000186092\"               |\n",
              "| \"ENST00000641515\"                     | \"ENSG00000186092\"               |\n",
              "+---------------------------------------+---------------------------------+\n",
              "\n",
              "+-------------------------------------+---------------------------------+\n",
              "| most_severe_consequence.gene_symbol | most_severe_consequence.biotype |\n",
              "+-------------------------------------+---------------------------------+\n",
              "| str                                 | str                             |\n",
              "+-------------------------------------+---------------------------------+\n",
              "| \"OR4F5\"                             | \"protein_coding\"                |\n",
              "| \"OR4F5\"                             | \"protein_coding\"                |\n",
              "| \"OR4F5\"                             | \"protein_coding\"                |\n",
              "| \"OR4F5\"                             | \"protein_coding\"                |\n",
              "| \"OR4F5\"                             | \"protein_coding\"                |\n",
              "| \"OR4F5\"                             | \"protein_coding\"                |\n",
              "| \"OR4F5\"                             | \"protein_coding\"                |\n",
              "| \"OR4F5\"                             | \"protein_coding\"                |\n",
              "| \"OR4F5\"                             | \"protein_coding\"                |\n",
              "| \"OR4F5\"                             | \"protein_coding\"                |\n",
              "+-------------------------------------+---------------------------------+\n",
              "\n",
              "+-------------------------------------------------+\n",
              "| most_severe_consequence.most_severe_consequence |\n",
              "+-------------------------------------------------+\n",
              "| str                                             |\n",
              "+-------------------------------------------------+\n",
              "| \"missense_variant\"                              |\n",
              "| \"missense_variant\"                              |\n",
              "| \"missense_variant\"                              |\n",
              "| \"missense_variant\"                              |\n",
              "| \"missense_variant\"                              |\n",
              "| \"missense_variant\"                              |\n",
              "| \"missense_variant\"                              |\n",
              "| \"missense_variant\"                              |\n",
              "| \"missense_variant\"                              |\n",
              "| \"missense_variant\"                              |\n",
              "+-------------------------------------------------+\n",
              "\n",
              "+-------------------------------------+-----------------------------------+\n",
              "| most_severe_consequence.mane_select | most_severe_consequence.canonical |\n",
              "+-------------------------------------+-----------------------------------+\n",
              "| str                                 |                             int32 |\n",
              "+-------------------------------------+-----------------------------------+\n",
              "| \"NM_001005484.2\"                    |                                 1 |\n",
              "| \"NM_001005484.2\"                    |                                 1 |\n",
              "| \"NM_001005484.2\"                    |                                 1 |\n",
              "| \"NM_001005484.2\"                    |                                 1 |\n",
              "| \"NM_001005484.2\"                    |                                 1 |\n",
              "| \"NM_001005484.2\"                    |                                 1 |\n",
              "| \"NM_001005484.2\"                    |                                 1 |\n",
              "| \"NM_001005484.2\"                    |                                 1 |\n",
              "| \"NM_001005484.2\"                    |                                 1 |\n",
              "| \"NM_001005484.2\"                    |                                 1 |\n",
              "+-------------------------------------+-----------------------------------+\n",
              "\n",
              "+-----------------------------+\n",
              "| most_severe_consequence.lof |\n",
              "+-----------------------------+\n",
              "| str                         |\n",
              "+-----------------------------+\n",
              "| NA                          |\n",
              "| NA                          |\n",
              "| NA                          |\n",
              "| NA                          |\n",
              "| NA                          |\n",
              "| NA                          |\n",
              "| NA                          |\n",
              "| NA                          |\n",
              "| NA                          |\n",
              "| NA                          |\n",
              "+-----------------------------+\n",
              "\n",
              "+-------------------------------------------------+\n",
              "| canonical_most_severe_consequence.transcript_id |\n",
              "+-------------------------------------------------+\n",
              "| str                                             |\n",
              "+-------------------------------------------------+\n",
              "| \"ENST00000641515\"                               |\n",
              "| \"ENST00000641515\"                               |\n",
              "| \"ENST00000641515\"                               |\n",
              "| \"ENST00000641515\"                               |\n",
              "| \"ENST00000641515\"                               |\n",
              "| \"ENST00000641515\"                               |\n",
              "| \"ENST00000641515\"                               |\n",
              "| \"ENST00000641515\"                               |\n",
              "| \"ENST00000641515\"                               |\n",
              "| \"ENST00000641515\"                               |\n",
              "+-------------------------------------------------+\n",
              "\n",
              "+-------------------------------------------+\n",
              "| canonical_most_severe_consequence.gene_id |\n",
              "+-------------------------------------------+\n",
              "| str                                       |\n",
              "+-------------------------------------------+\n",
              "| \"ENSG00000186092\"                         |\n",
              "| \"ENSG00000186092\"                         |\n",
              "| \"ENSG00000186092\"                         |\n",
              "| \"ENSG00000186092\"                         |\n",
              "| \"ENSG00000186092\"                         |\n",
              "| \"ENSG00000186092\"                         |\n",
              "| \"ENSG00000186092\"                         |\n",
              "| \"ENSG00000186092\"                         |\n",
              "| \"ENSG00000186092\"                         |\n",
              "| \"ENSG00000186092\"                         |\n",
              "+-------------------------------------------+\n",
              "\n",
              "+-----------------------------------------------+\n",
              "| canonical_most_severe_consequence.gene_symbol |\n",
              "+-----------------------------------------------+\n",
              "| str                                           |\n",
              "+-----------------------------------------------+\n",
              "| \"OR4F5\"                                       |\n",
              "| \"OR4F5\"                                       |\n",
              "| \"OR4F5\"                                       |\n",
              "| \"OR4F5\"                                       |\n",
              "| \"OR4F5\"                                       |\n",
              "| \"OR4F5\"                                       |\n",
              "| \"OR4F5\"                                       |\n",
              "| \"OR4F5\"                                       |\n",
              "| \"OR4F5\"                                       |\n",
              "| \"OR4F5\"                                       |\n",
              "+-----------------------------------------------+\n",
              "\n",
              "+-------------------------------------------+\n",
              "| canonical_most_severe_consequence.biotype |\n",
              "+-------------------------------------------+\n",
              "| str                                       |\n",
              "+-------------------------------------------+\n",
              "| \"protein_coding\"                          |\n",
              "| \"protein_coding\"                          |\n",
              "| \"protein_coding\"                          |\n",
              "| \"protein_coding\"                          |\n",
              "| \"protein_coding\"                          |\n",
              "| \"protein_coding\"                          |\n",
              "| \"protein_coding\"                          |\n",
              "| \"protein_coding\"                          |\n",
              "| \"protein_coding\"                          |\n",
              "| \"protein_coding\"                          |\n",
              "+-------------------------------------------+\n",
              "\n",
              "+-----------------------------------------------------------+\n",
              "| canonical_most_severe_consequence.most_severe_consequence |\n",
              "+-----------------------------------------------------------+\n",
              "| str                                                       |\n",
              "+-----------------------------------------------------------+\n",
              "| \"missense_variant\"                                        |\n",
              "| \"missense_variant\"                                        |\n",
              "| \"missense_variant\"                                        |\n",
              "| \"missense_variant\"                                        |\n",
              "| \"missense_variant\"                                        |\n",
              "| \"missense_variant\"                                        |\n",
              "| \"missense_variant\"                                        |\n",
              "| \"missense_variant\"                                        |\n",
              "| \"missense_variant\"                                        |\n",
              "| \"missense_variant\"                                        |\n",
              "+-----------------------------------------------------------+\n",
              "\n",
              "+-----------------------------------------------+\n",
              "| canonical_most_severe_consequence.mane_select |\n",
              "+-----------------------------------------------+\n",
              "| str                                           |\n",
              "+-----------------------------------------------+\n",
              "| \"NM_001005484.2\"                              |\n",
              "| \"NM_001005484.2\"                              |\n",
              "| \"NM_001005484.2\"                              |\n",
              "| \"NM_001005484.2\"                              |\n",
              "| \"NM_001005484.2\"                              |\n",
              "| \"NM_001005484.2\"                              |\n",
              "| \"NM_001005484.2\"                              |\n",
              "| \"NM_001005484.2\"                              |\n",
              "| \"NM_001005484.2\"                              |\n",
              "| \"NM_001005484.2\"                              |\n",
              "+-----------------------------------------------+\n",
              "\n",
              "+---------------------------------------------+\n",
              "| canonical_most_severe_consequence.canonical |\n",
              "+---------------------------------------------+\n",
              "|                                       int32 |\n",
              "+---------------------------------------------+\n",
              "|                                           1 |\n",
              "|                                           1 |\n",
              "|                                           1 |\n",
              "|                                           1 |\n",
              "|                                           1 |\n",
              "|                                           1 |\n",
              "|                                           1 |\n",
              "|                                           1 |\n",
              "|                                           1 |\n",
              "|                                           1 |\n",
              "+---------------------------------------------+\n",
              "\n",
              "+---------------------------------------+\n",
              "| canonical_most_severe_consequence.lof |\n",
              "+---------------------------------------+\n",
              "| str                                   |\n",
              "+---------------------------------------+\n",
              "| NA                                    |\n",
              "| NA                                    |\n",
              "| NA                                    |\n",
              "| NA                                    |\n",
              "| NA                                    |\n",
              "| NA                                    |\n",
              "| NA                                    |\n",
              "| NA                                    |\n",
              "| NA                                    |\n",
              "| NA                                    |\n",
              "+---------------------------------------+\n",
              "\n",
              "+---------------------------------------------------+\n",
              "| mane_select_most_severe_consequence.transcript_id |\n",
              "+---------------------------------------------------+\n",
              "| str                                               |\n",
              "+---------------------------------------------------+\n",
              "| \"ENST00000641515\"                                 |\n",
              "| \"ENST00000641515\"                                 |\n",
              "| \"ENST00000641515\"                                 |\n",
              "| \"ENST00000641515\"                                 |\n",
              "| \"ENST00000641515\"                                 |\n",
              "| \"ENST00000641515\"                                 |\n",
              "| \"ENST00000641515\"                                 |\n",
              "| \"ENST00000641515\"                                 |\n",
              "| \"ENST00000641515\"                                 |\n",
              "| \"ENST00000641515\"                                 |\n",
              "+---------------------------------------------------+\n",
              "\n",
              "+---------------------------------------------+\n",
              "| mane_select_most_severe_consequence.gene_id |\n",
              "+---------------------------------------------+\n",
              "| str                                         |\n",
              "+---------------------------------------------+\n",
              "| \"ENSG00000186092\"                           |\n",
              "| \"ENSG00000186092\"                           |\n",
              "| \"ENSG00000186092\"                           |\n",
              "| \"ENSG00000186092\"                           |\n",
              "| \"ENSG00000186092\"                           |\n",
              "| \"ENSG00000186092\"                           |\n",
              "| \"ENSG00000186092\"                           |\n",
              "| \"ENSG00000186092\"                           |\n",
              "| \"ENSG00000186092\"                           |\n",
              "| \"ENSG00000186092\"                           |\n",
              "+---------------------------------------------+\n",
              "\n",
              "+-------------------------------------------------+\n",
              "| mane_select_most_severe_consequence.gene_symbol |\n",
              "+-------------------------------------------------+\n",
              "| str                                             |\n",
              "+-------------------------------------------------+\n",
              "| \"OR4F5\"                                         |\n",
              "| \"OR4F5\"                                         |\n",
              "| \"OR4F5\"                                         |\n",
              "| \"OR4F5\"                                         |\n",
              "| \"OR4F5\"                                         |\n",
              "| \"OR4F5\"                                         |\n",
              "| \"OR4F5\"                                         |\n",
              "| \"OR4F5\"                                         |\n",
              "| \"OR4F5\"                                         |\n",
              "| \"OR4F5\"                                         |\n",
              "+-------------------------------------------------+\n",
              "\n",
              "+---------------------------------------------+\n",
              "| mane_select_most_severe_consequence.biotype |\n",
              "+---------------------------------------------+\n",
              "| str                                         |\n",
              "+---------------------------------------------+\n",
              "| \"protein_coding\"                            |\n",
              "| \"protein_coding\"                            |\n",
              "| \"protein_coding\"                            |\n",
              "| \"protein_coding\"                            |\n",
              "| \"protein_coding\"                            |\n",
              "| \"protein_coding\"                            |\n",
              "| \"protein_coding\"                            |\n",
              "| \"protein_coding\"                            |\n",
              "| \"protein_coding\"                            |\n",
              "| \"protein_coding\"                            |\n",
              "+---------------------------------------------+\n",
              "\n",
              "+-------------------------------------------------------------+\n",
              "| mane_select_most_severe_consequence.most_severe_consequence |\n",
              "+-------------------------------------------------------------+\n",
              "| str                                                         |\n",
              "+-------------------------------------------------------------+\n",
              "| \"missense_variant\"                                          |\n",
              "| \"missense_variant\"                                          |\n",
              "| \"missense_variant\"                                          |\n",
              "| \"missense_variant\"                                          |\n",
              "| \"missense_variant\"                                          |\n",
              "| \"missense_variant\"                                          |\n",
              "| \"missense_variant\"                                          |\n",
              "| \"missense_variant\"                                          |\n",
              "| \"missense_variant\"                                          |\n",
              "| \"missense_variant\"                                          |\n",
              "+-------------------------------------------------------------+\n",
              "\n",
              "+-------------------------------------------------+\n",
              "| mane_select_most_severe_consequence.mane_select |\n",
              "+-------------------------------------------------+\n",
              "| str                                             |\n",
              "+-------------------------------------------------+\n",
              "| \"NM_001005484.2\"                                |\n",
              "| \"NM_001005484.2\"                                |\n",
              "| \"NM_001005484.2\"                                |\n",
              "| \"NM_001005484.2\"                                |\n",
              "| \"NM_001005484.2\"                                |\n",
              "| \"NM_001005484.2\"                                |\n",
              "| \"NM_001005484.2\"                                |\n",
              "| \"NM_001005484.2\"                                |\n",
              "| \"NM_001005484.2\"                                |\n",
              "| \"NM_001005484.2\"                                |\n",
              "+-------------------------------------------------+\n",
              "\n",
              "+-----------------------------------------------+\n",
              "| mane_select_most_severe_consequence.canonical |\n",
              "+-----------------------------------------------+\n",
              "|                                         int32 |\n",
              "+-----------------------------------------------+\n",
              "|                                             1 |\n",
              "|                                             1 |\n",
              "|                                             1 |\n",
              "|                                             1 |\n",
              "|                                             1 |\n",
              "|                                             1 |\n",
              "|                                             1 |\n",
              "|                                             1 |\n",
              "|                                             1 |\n",
              "|                                             1 |\n",
              "+-----------------------------------------------+\n",
              "\n",
              "+-----------------------------------------+-------+----------+----------+\n",
              "| mane_select_most_severe_consequence.lof |    AC |       AF | ac_group |\n",
              "+-----------------------------------------+-------+----------+----------+\n",
              "| str                                     | int32 |  float64 | str      |\n",
              "+-----------------------------------------+-------+----------+----------+\n",
              "| NA                                      |     0 | 0.00e+00 | \"0\"      |\n",
              "| NA                                      |     0 | 0.00e+00 | \"0\"      |\n",
              "| NA                                      |     0 | 0.00e+00 | \"0\"      |\n",
              "| NA                                      |     0 | 0.00e+00 | \"0\"      |\n",
              "| NA                                      |     0 | 0.00e+00 | \"0\"      |\n",
              "| NA                                      |     0 | 0.00e+00 | \"0\"      |\n",
              "| NA                                      |     0 | 0.00e+00 | \"0\"      |\n",
              "| NA                                      |     0 | 0.00e+00 | \"0\"      |\n",
              "| NA                                      |     0 | 0.00e+00 | \"0\"      |\n",
              "| NA                                      |     0 | 0.00e+00 | \"0\"      |\n",
              "+-----------------------------------------+-------+----------+----------+\n",
              "\n",
              "+--------------+\n",
              "| downsampling |\n",
              "+--------------+\n",
              "| str          |\n",
              "+--------------+\n",
              "| \"AC_ds10\"    |\n",
              "| \"AC_ds100\"   |\n",
              "| \"AC_ds500\"   |\n",
              "| \"AC_ds1000\"  |\n",
              "| \"AC_ds2000\"  |\n",
              "| \"AC_ds2884\"  |\n",
              "| \"AC_ds5000\"  |\n",
              "| \"AC_ds10000\" |\n",
              "| \"AC_ds13068\" |\n",
              "| \"AC_ds16740\" |\n",
              "+--------------+\n",
              "showing top 10 rows"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        }
      ],
      "source": [
        "from typing import Optional, Union\n",
        "\n",
        "from gnomad.utils.vep import CSQ_ORDER, LOF_CSQ_SET, get_most_severe_consequence_expr\n",
        "\n",
        "def create_ac_group_expr(ac_expr):\n",
        "    return (\n",
        "        hl.case()\n",
        "        .when(ac_expr == 0, \"0\")\n",
        "        .when(ac_expr == 1, \"1\") \n",
        "        .when(ac_expr == 2, \"2\") \n",
        "        .when(ac_expr == 3, \"3\") \n",
        "        .when(ac_expr == 4, \"4\") \n",
        "        .when(ac_expr == 5, \"5\") \n",
        "        .when(ac_expr == 6, \"6\") \n",
        "        .when(ac_expr == 7, \"7\") \n",
        "        .when(ac_expr == 8, \"8\") \n",
        "        .when(ac_expr == 9, \"9\") \n",
        "        .when(ac_expr == 10, \"10\") \n",
        "        .when((ac_expr > 10) & (ac_expr <= 100), \"10-100\") \n",
        "        .when((ac_expr > 100) & (ac_expr <= 1000), \"100-1000\") \n",
        "        .when((ac_expr > 1000) & (ac_expr <= 10000), \"1000-10000\") \n",
        "        .when(ac_expr > 10000, \">10000\") \n",
        "        .or_missing()\n",
        "    )\n",
        "    \n",
        "\n",
        "def get_most_severe_consequence_expr(\n",
        "    csq_expr: hl.expr.ArrayExpression,\n",
        "    csq_order: Optional[List[str]] = None,\n",
        "    csq_field: str = \"most_severe_consequence\",\n",
        ") -> Union[hl.expr.StringExpression, hl.expr.StructExpression]:\n",
        "    \"\"\"\n",
        "    Get the most severe consequence from a collection of consequences.\n",
        "\n",
        "    This is for a given transcript, as there are often multiple annotations for a single\n",
        "    transcript: e.g. splice_region_variant&intron_variant -> splice_region_variant\n",
        "\n",
        "    :param csq_expr: ArrayExpression of consequences.\n",
        "    :param csq_order: Optional list indicating the order of VEP consequences, sorted\n",
        "        from high to low impact. Default is None, which uses the value of the\n",
        "        `CSQ_ORDER` global.\n",
        "    :return: Most severe consequence in `csq_expr`.\n",
        "    \"\"\"\n",
        "    if csq_order is None:\n",
        "        csq_order = CSQ_ORDER\n",
        "    csqs = hl.literal(csq_order)\n",
        "\n",
        "    \n",
        "\n",
        "    if csq_expr.dtype == hl.tarray(hl.tstr):\n",
        "        return csqs.find(lambda c: csq_expr.contains(c))\n",
        "    else:\n",
        "        ms_csq = csqs.find(lambda c: csq_expr.map(lambda x: x[csq_field]).contains(c))\n",
        "        return csq_expr.filter(lambda x: x[csq_field] == ms_csq).first()\n",
        "\n",
        "\n",
        "\n",
        "downsampling_ht_path = \"gs://gnomad/v4.1/constraint_coverage_corrected/preprocessed_data/gnomad.v4.1.context.preprocessed.ht\"\n",
        "\n",
        "csq_to_keep = [\"missense_variant\", \"synonymous_variant\"] + list(LOF_CSQ_SET)\n",
        "\n",
        "ht = hl.read_table(downsampling_ht_path)\n",
        "freq_meta = hl.eval(ht.exomes_freq_meta)\n",
        "ht = ht.filter(\n",
        "    hl.is_defined(ht.filters.exomes) \n",
        "    & (ht.filters.exomes.length() == 0)\n",
        "    & (hl.set(csq_to_keep).contains(ht.vep.most_severe_consequence))\n",
        ")\n",
        "\n",
        "ht = ht.select(\n",
        "    #\"context\",\n",
        "    #'transition',\n",
        "    #'cpg',\n",
        "    'mutation_type',\n",
        "    #'mutation_type_model',\n",
        "    #'exomes_coverage',\n",
        "    transcript_consequences=ht.vep.transcript_consequences.map(\n",
        "        lambda x: x.select(\n",
        "            \"transcript_id\", \n",
        "            \"gene_id\", \n",
        "            \"gene_symbol\", \n",
        "            \"biotype\", \n",
        "            \"most_severe_consequence\", \n",
        "            \"mane_select\", \n",
        "            \"canonical\", \n",
        "            \"lof\", \n",
        "        )\n",
        "    ).filter(\n",
        "        lambda x: (\n",
        "            x.transcript_id.startswith(\"ENST\") \n",
        "            & (x.biotype == \"protein_coding\")\n",
        "            & hl.set(csq_to_keep).contains(x.most_severe_consequence)\n",
        "        )\n",
        "    ),\n",
        "    freq=ht.calibrate_mu.exomes_freq.map(\n",
        "        lambda x: x.select(\n",
        "            \"AC\", \n",
        "            \"AF\",\n",
        "            ac_group=create_ac_group_expr(x.AC)\n",
        "        )\n",
        "    )\n",
        ")\n",
        "ht = ht.filter(ht.transcript_consequences.length() > 0)\n",
        "ht = ht.select(\n",
        "    'mutation_type',\n",
        "    freq=[\n",
        "        ht.freq[i].annotate(downsampling=f\"AC_ds{m['downsampling']}\") \n",
        "        for i, m in enumerate(freq_meta)\n",
        "        if \"downsampling\" in m.keys() and \"gen_anc\" in m.keys() and m[\"gen_anc\"] == \"global\"\n",
        "    ],\n",
        "    most_severe_consequence=get_most_severe_consequence_expr(\n",
        "        ht.transcript_consequences\n",
        "    ),\n",
        "    canonical_most_severe_consequence=get_most_severe_consequence_expr(\n",
        "        ht.transcript_consequences.filter(lambda x: hl.or_else(x.canonical==1, False))\n",
        "    ),\n",
        "    mane_select_most_severe_consequence=get_most_severe_consequence_expr(\n",
        "        ht.transcript_consequences.filter(lambda x: hl.is_defined(x.mane_select))\n",
        "    )\n",
        ")\n",
        "ht = ht.explode(ht.freq)\n",
        "ht = ht.annotate(**ht.freq).drop(\"freq\")\n",
        "ht = ht.checkpoint(\n",
        "    \"gs://gnomad-tmp-4day/julia/gnomad.v4.1.constraint.preprocessed.for_ac_group_computation.ht\",\n",
        "    #_read_if_exists=True,\n",
        "    overwrite=True\n",
        ")\n",
        "\n",
        "print(\"\\nTable schema:\")\n",
        "ht.describe()\n",
        "ht.show()\n"
      ]
    },
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {},
      "outputs": [
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\"></pre>\n"
            ],
            "text/plain": [
              "\u001b[?25l"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "application/vnd.jupyter.widget-view+json": {
              "model_id": "cd906d038fcb4d73b61ef38218bc04a4",
              "version_major": 2,
              "version_minor": 0
            },
            "text/plain": [
              "Output()"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\">\n",
              "</pre>\n"
            ],
            "text/plain": [
              "\n",
              "\u001b[?25h"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\"></pre>\n"
            ],
            "text/plain": [
              "\u001b[?25l"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "application/vnd.jupyter.widget-view+json": {
              "model_id": "c66fa2a22f51424bb467ae7f1ac8624f",
              "version_major": 2,
              "version_minor": 0
            },
            "text/plain": [
              "Output()"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\">\n",
              "</pre>\n"
            ],
            "text/plain": [
              "\n",
              "\u001b[?25h"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\"></pre>\n"
            ],
            "text/plain": [
              "\u001b[?25l"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "application/vnd.jupyter.widget-view+json": {
              "model_id": "2300c14443754ff89e5c01e8243b9d28",
              "version_major": 2,
              "version_minor": 0
            },
            "text/plain": [
              "Output()"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "name": "stderr",
          "output_type": "stream",
          "text": [
            "WARNING (hailtop.utils 828): A transient error occured. We will automatically retry. Do not be alarmed. We have thus far seen 2 transient errors (next delay: 3.653s). The most recent error was <class 'aiohttp.client_exceptions.ClientConnectorError'> Cannot connect to host batch.hail.is:443 ssl:default [Connect call failed ('35.188.91.25', 443)]. \n"
          ]
        },
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\">\r\u001b[2KReceived a keyboard interrupt, cancelling the batch...\n",
              "execute(...): <a href=\"https://batch.hail.is/batches/8325882\" target=\"_blank\">8325882</a> <span style=\"color: #729c1f; text-decoration-color: #729c1f\">━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╸</span> <span style=\"color: #800080; text-decoration-color: #800080\">100%</span> <span style=\"color: #008000; text-decoration-color: #008000\">68677/68679</span> <span style=\"color: #008080; text-decoration-color: #008080\">-:--:--</span> <span style=\"color: #808000; text-decoration-color: #808000\">5:39:25</span></pre>\n"
            ],
            "text/plain": [
              "\r\u001b[2KReceived a keyboard interrupt, cancelling the batch...\n",
              "execute(...): \u001b]8;id=104996;https://batch.hail.is/batches/8325882\u001b\\8325882\u001b]8;;\u001b\\ \u001b[38;2;114;156;31m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\u001b[0m\u001b[38;2;114;156;31m╸\u001b[0m \u001b[35m100%\u001b[0m \u001b[32m68677/68679\u001b[0m \u001b[36m-:--:--\u001b[0m \u001b[33m5:39:25\u001b[0m"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "data": {
            "text/html": [
              "<pre style=\"white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace\">\n",
              "</pre>\n"
            ],
            "text/plain": [
              "\n",
              "\u001b[?25h"
            ]
          },
          "metadata": {},
          "output_type": "display_data"
        },
        {
          "ename": "KeyboardInterrupt",
          "evalue": "",
          "output_type": "error",
          "traceback": [
            "\u001b[31m---------------------------------------------------------------------------\u001b[39m",
            "\u001b[31mKeyboardInterrupt\u001b[39m                         Traceback (most recent call last)",
            "\u001b[36mCell\u001b[39m\u001b[36m \u001b[39m\u001b[32mIn[16]\u001b[39m\u001b[32m, line 18\u001b[39m\n\u001b[32m      1\u001b[39m ht = hl.read_table(\u001b[33m\"\u001b[39m\u001b[33mgs://gnomad-tmp-4day/julia/gnomad.v4.1.constraint.preprocessed.for_ac_group_computation.ht\u001b[39m\u001b[33m\"\u001b[39m)\n\u001b[32m      2\u001b[39m agg_ht = \u001b[43mht\u001b[49m\u001b[43m.\u001b[49m\u001b[43mgroup_by\u001b[49m\u001b[43m(\u001b[49m\u001b[33;43m\"\u001b[39;49m\u001b[33;43mmutation_type\u001b[39;49m\u001b[33;43m\"\u001b[39;49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[33;43m\"\u001b[39;49m\u001b[33;43mdownsampling\u001b[39;49m\u001b[33;43m\"\u001b[39;49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[33;43m\"\u001b[39;49m\u001b[33;43mac_group\u001b[39;49m\u001b[33;43m\"\u001b[39;49m\u001b[43m)\u001b[49m\u001b[43m.\u001b[49m\u001b[43maggregate\u001b[49m\u001b[43m(\u001b[49m\n\u001b[32m      3\u001b[39m \u001b[43m    \u001b[49m\u001b[43mn\u001b[49m\u001b[43m=\u001b[49m\u001b[43mhl\u001b[49m\u001b[43m.\u001b[49m\u001b[43mstruct\u001b[49m\u001b[43m(\u001b[49m\n\u001b[32m      4\u001b[39m \u001b[43m        \u001b[49m\u001b[43mtotal\u001b[49m\u001b[43m=\u001b[49m\u001b[43mhl\u001b[49m\u001b[43m.\u001b[49m\u001b[43magg\u001b[49m\u001b[43m.\u001b[49m\u001b[43mcount\u001b[49m\u001b[43m(\u001b[49m\u001b[43m)\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m      5\u001b[39m \u001b[43m        \u001b[49m\u001b[43mmost_severe_consequence\u001b[49m\u001b[43m=\u001b[49m\u001b[43mhl\u001b[49m\u001b[43m.\u001b[49m\u001b[43magg\u001b[49m\u001b[43m.\u001b[49m\u001b[43mgroup_by\u001b[49m\u001b[43m(\u001b[49m\n\u001b[32m      6\u001b[39m \u001b[43m            \u001b[49m\u001b[43mht\u001b[49m\u001b[43m.\u001b[49m\u001b[43mmost_severe_consequence\u001b[49m\u001b[43m.\u001b[49m\u001b[43mselect\u001b[49m\u001b[43m(\u001b[49m\u001b[33;43m\"\u001b[39;49m\u001b[33;43mmost_severe_consequence\u001b[39;49m\u001b[33;43m\"\u001b[39;49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[33;43m\"\u001b[39;49m\u001b[33;43mlof\u001b[39;49m\u001b[33;43m\"\u001b[39;49m\u001b[43m)\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m      7\u001b[39m \u001b[43m            \u001b[49m\u001b[43mhl\u001b[49m\u001b[43m.\u001b[49m\u001b[43magg\u001b[49m\u001b[43m.\u001b[49m\u001b[43mcount\u001b[49m\u001b[43m(\u001b[49m\u001b[43m)\u001b[49m\n\u001b[32m      8\u001b[39m \u001b[43m        \u001b[49m\u001b[43m)\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m      9\u001b[39m \u001b[43m        \u001b[49m\u001b[43mcanonical_most_severe_consequence\u001b[49m\u001b[43m=\u001b[49m\u001b[43mhl\u001b[49m\u001b[43m.\u001b[49m\u001b[43magg\u001b[49m\u001b[43m.\u001b[49m\u001b[43mgroup_by\u001b[49m\u001b[43m(\u001b[49m\n\u001b[32m     10\u001b[39m \u001b[43m            \u001b[49m\u001b[43mht\u001b[49m\u001b[43m.\u001b[49m\u001b[43mcanonical_most_severe_consequence\u001b[49m\u001b[43m.\u001b[49m\u001b[43mselect\u001b[49m\u001b[43m(\u001b[49m\u001b[33;43m\"\u001b[39;49m\u001b[33;43mmost_severe_consequence\u001b[39;49m\u001b[33;43m\"\u001b[39;49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[33;43m\"\u001b[39;49m\u001b[33;43mlof\u001b[39;49m\u001b[33;43m\"\u001b[39;49m\u001b[43m)\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m     11\u001b[39m \u001b[43m            \u001b[49m\u001b[43mhl\u001b[49m\u001b[43m.\u001b[49m\u001b[43magg\u001b[49m\u001b[43m.\u001b[49m\u001b[43mcount\u001b[49m\u001b[43m(\u001b[49m\u001b[43m)\u001b[49m\n\u001b[32m     12\u001b[39m \u001b[43m        \u001b[49m\u001b[43m)\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m     13\u001b[39m \u001b[43m        \u001b[49m\u001b[43mmane_select_most_severe_consequence\u001b[49m\u001b[43m=\u001b[49m\u001b[43mhl\u001b[49m\u001b[43m.\u001b[49m\u001b[43magg\u001b[49m\u001b[43m.\u001b[49m\u001b[43mgroup_by\u001b[49m\u001b[43m(\u001b[49m\n\u001b[32m     14\u001b[39m \u001b[43m            \u001b[49m\u001b[43mht\u001b[49m\u001b[43m.\u001b[49m\u001b[43mmane_select_most_severe_consequence\u001b[49m\u001b[43m.\u001b[49m\u001b[43mselect\u001b[49m\u001b[43m(\u001b[49m\u001b[33;43m\"\u001b[39;49m\u001b[33;43mmost_severe_consequence\u001b[39;49m\u001b[33;43m\"\u001b[39;49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[33;43m\"\u001b[39;49m\u001b[33;43mlof\u001b[39;49m\u001b[33;43m\"\u001b[39;49m\u001b[43m)\u001b[49m\u001b[43m,\u001b[49m\n\u001b[32m     15\u001b[39m \u001b[43m            \u001b[49m\u001b[43mhl\u001b[49m\u001b[43m.\u001b[49m\u001b[43magg\u001b[49m\u001b[43m.\u001b[49m\u001b[43mcount\u001b[49m\u001b[43m(\u001b[49m\u001b[43m)\u001b[49m\n\u001b[32m     16\u001b[39m \u001b[43m        \u001b[49m\u001b[43m)\u001b[49m\n\u001b[32m     17\u001b[39m \u001b[43m    \u001b[49m\u001b[43m)\u001b[49m\n\u001b[32m---> \u001b[39m\u001b[32m18\u001b[39m \u001b[43m)\u001b[49m\u001b[43m.\u001b[49m\u001b[43mcheckpoint\u001b[49m\u001b[43m(\u001b[49m\n\u001b[32m     19\u001b[39m \u001b[43m    \u001b[49m\u001b[33;43m\"\u001b[39;49m\u001b[33;43mgs://gnomad-tmp-4day/julia/gnomad.v4.1.constraint.preprocessed.for_ac_group_computation.agg_ht.ht\u001b[39;49m\u001b[33;43m\"\u001b[39;49m\u001b[43m,\u001b[49m\n\u001b[32m     20\u001b[39m \u001b[43m    \u001b[49m\u001b[43m_read_if_exists\u001b[49m\u001b[43m=\u001b[49m\u001b[38;5;28;43;01mTrue\u001b[39;49;00m\u001b[43m,\u001b[49m\n\u001b[32m     21\u001b[39m \u001b[43m    \u001b[49m\u001b[38;5;66;43;03m#overwrite=True\u001b[39;49;00m\n\u001b[32m     22\u001b[39m \u001b[43m)\u001b[49m\n\u001b[32m     23\u001b[39m agg_ht.show()\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/decorator.py:235\u001b[39m, in \u001b[36mdecorate.<locals>.fun\u001b[39m\u001b[34m(*args, **kw)\u001b[39m\n\u001b[32m    233\u001b[39m \u001b[38;5;28;01mif\u001b[39;00m \u001b[38;5;129;01mnot\u001b[39;00m kwsyntax:\n\u001b[32m    234\u001b[39m     args, kw = fix(args, kw, sig)\n\u001b[32m--> \u001b[39m\u001b[32m235\u001b[39m \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[43mcaller\u001b[49m\u001b[43m(\u001b[49m\u001b[43mfunc\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43m*\u001b[49m\u001b[43m(\u001b[49m\u001b[43mextras\u001b[49m\u001b[43m \u001b[49m\u001b[43m+\u001b[49m\u001b[43m \u001b[49m\u001b[43margs\u001b[49m\u001b[43m)\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43m*\u001b[49m\u001b[43m*\u001b[49m\u001b[43mkw\u001b[49m\u001b[43m)\u001b[49m\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/hail/typecheck/check.py:585\u001b[39m, in \u001b[36m_make_dec.<locals>.wrapper\u001b[39m\u001b[34m(__original_func, *args, **kwargs)\u001b[39m\n\u001b[32m    582\u001b[39m \u001b[38;5;129m@decorator\u001b[39m\n\u001b[32m    583\u001b[39m \u001b[38;5;28;01mdef\u001b[39;00m\u001b[38;5;250m \u001b[39m\u001b[34mwrapper\u001b[39m(__original_func: Callable[..., T], *args, **kwargs) -> T:\n\u001b[32m    584\u001b[39m     args_, kwargs_ = check_all(__original_func, args, kwargs, checkers, is_method=is_method)\n\u001b[32m--> \u001b[39m\u001b[32m585\u001b[39m     \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[43m__original_func\u001b[49m\u001b[43m(\u001b[49m\u001b[43m*\u001b[49m\u001b[43margs_\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43m*\u001b[49m\u001b[43m*\u001b[49m\u001b[43mkwargs_\u001b[49m\u001b[43m)\u001b[49m\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/hail/table.py:1963\u001b[39m, in \u001b[36mTable.checkpoint\u001b[39m\u001b[34m(self, output, overwrite, stage_locally, _codec_spec, _read_if_exists, _intervals, _filter_intervals)\u001b[39m\n\u001b[32m   1960\u001b[39m hl.current_backend().validate_file(output)\n\u001b[32m   1962\u001b[39m \u001b[38;5;28;01mif\u001b[39;00m \u001b[38;5;129;01mnot\u001b[39;00m _read_if_exists \u001b[38;5;129;01mor\u001b[39;00m \u001b[38;5;129;01mnot\u001b[39;00m hl.hadoop_exists(\u001b[33mf\u001b[39m\u001b[33m'\u001b[39m\u001b[38;5;132;01m{\u001b[39;00moutput\u001b[38;5;132;01m}\u001b[39;00m\u001b[33m/_SUCCESS\u001b[39m\u001b[33m'\u001b[39m):\n\u001b[32m-> \u001b[39m\u001b[32m1963\u001b[39m     \u001b[38;5;28;43mself\u001b[39;49m\u001b[43m.\u001b[49m\u001b[43mwrite\u001b[49m\u001b[43m(\u001b[49m\u001b[43moutput\u001b[49m\u001b[43m=\u001b[49m\u001b[43moutput\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43moverwrite\u001b[49m\u001b[43m=\u001b[49m\u001b[43moverwrite\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mstage_locally\u001b[49m\u001b[43m=\u001b[49m\u001b[43mstage_locally\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43m_codec_spec\u001b[49m\u001b[43m=\u001b[49m\u001b[43m_codec_spec\u001b[49m\u001b[43m)\u001b[49m\n\u001b[32m   1964\u001b[39m     _assert_type = \u001b[38;5;28mself\u001b[39m._type\n\u001b[32m   1965\u001b[39m     _load_refs = \u001b[38;5;28;01mFalse\u001b[39;00m\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/decorator.py:235\u001b[39m, in \u001b[36mdecorate.<locals>.fun\u001b[39m\u001b[34m(*args, **kw)\u001b[39m\n\u001b[32m    233\u001b[39m \u001b[38;5;28;01mif\u001b[39;00m \u001b[38;5;129;01mnot\u001b[39;00m kwsyntax:\n\u001b[32m    234\u001b[39m     args, kw = fix(args, kw, sig)\n\u001b[32m--> \u001b[39m\u001b[32m235\u001b[39m \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[43mcaller\u001b[49m\u001b[43m(\u001b[49m\u001b[43mfunc\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43m*\u001b[49m\u001b[43m(\u001b[49m\u001b[43mextras\u001b[49m\u001b[43m \u001b[49m\u001b[43m+\u001b[49m\u001b[43m \u001b[49m\u001b[43margs\u001b[49m\u001b[43m)\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43m*\u001b[49m\u001b[43m*\u001b[49m\u001b[43mkw\u001b[49m\u001b[43m)\u001b[49m\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/hail/typecheck/check.py:585\u001b[39m, in \u001b[36m_make_dec.<locals>.wrapper\u001b[39m\u001b[34m(__original_func, *args, **kwargs)\u001b[39m\n\u001b[32m    582\u001b[39m \u001b[38;5;129m@decorator\u001b[39m\n\u001b[32m    583\u001b[39m \u001b[38;5;28;01mdef\u001b[39;00m\u001b[38;5;250m \u001b[39m\u001b[34mwrapper\u001b[39m(__original_func: Callable[..., T], *args, **kwargs) -> T:\n\u001b[32m    584\u001b[39m     args_, kwargs_ = check_all(__original_func, args, kwargs, checkers, is_method=is_method)\n\u001b[32m--> \u001b[39m\u001b[32m585\u001b[39m     \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[43m__original_func\u001b[49m\u001b[43m(\u001b[49m\u001b[43m*\u001b[49m\u001b[43margs_\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43m*\u001b[49m\u001b[43m*\u001b[49m\u001b[43mkwargs_\u001b[49m\u001b[43m)\u001b[49m\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/hail/table.py:2005\u001b[39m, in \u001b[36mTable.write\u001b[39m\u001b[34m(self, output, overwrite, stage_locally, _codec_spec)\u001b[39m\n\u001b[32m   1979\u001b[39m \u001b[38;5;250m\u001b[39m\u001b[33;03m\"\"\"Write to disk.\u001b[39;00m\n\u001b[32m   1980\u001b[39m \n\u001b[32m   1981\u001b[39m \u001b[33;03mExamples\u001b[39;00m\n\u001b[32m   (...)\u001b[39m\u001b[32m   2000\u001b[39m \u001b[33;03m    If ``True``, overwrite an existing file at the destination.\u001b[39;00m\n\u001b[32m   2001\u001b[39m \u001b[33;03m\"\"\"\u001b[39;00m\n\u001b[32m   2003\u001b[39m hl.current_backend().validate_file(output)\n\u001b[32m-> \u001b[39m\u001b[32m2005\u001b[39m \u001b[43mEnv\u001b[49m\u001b[43m.\u001b[49m\u001b[43mbackend\u001b[49m\u001b[43m(\u001b[49m\u001b[43m)\u001b[49m\u001b[43m.\u001b[49m\u001b[43mexecute\u001b[49m\u001b[43m(\u001b[49m\n\u001b[32m   2006\u001b[39m \u001b[43m    \u001b[49m\u001b[43mir\u001b[49m\u001b[43m.\u001b[49m\u001b[43mTableWrite\u001b[49m\u001b[43m(\u001b[49m\u001b[38;5;28;43mself\u001b[39;49m\u001b[43m.\u001b[49m\u001b[43m_tir\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mir\u001b[49m\u001b[43m.\u001b[49m\u001b[43mTableNativeWriter\u001b[49m\u001b[43m(\u001b[49m\u001b[43moutput\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43moverwrite\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mstage_locally\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43m_codec_spec\u001b[49m\u001b[43m)\u001b[49m\u001b[43m)\u001b[49m\n\u001b[32m   2007\u001b[39m \u001b[43m\u001b[49m\u001b[43m)\u001b[49m\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/hail/backend/backend.py:179\u001b[39m, in \u001b[36mBackend.execute\u001b[39m\u001b[34m(self, ir, timed)\u001b[39m\n\u001b[32m    177\u001b[39m payload = ExecutePayload(\u001b[38;5;28mself\u001b[39m._render_ir(ir), \u001b[33m'\u001b[39m\u001b[33m{\u001b[39m\u001b[33m\"\u001b[39m\u001b[33mname\u001b[39m\u001b[33m\"\u001b[39m\u001b[33m:\u001b[39m\u001b[33m\"\u001b[39m\u001b[33mStreamBufferSpec\u001b[39m\u001b[33m\"\u001b[39m\u001b[33m}\u001b[39m\u001b[33m'\u001b[39m, timed)\n\u001b[32m    178\u001b[39m \u001b[38;5;28;01mtry\u001b[39;00m:\n\u001b[32m--> \u001b[39m\u001b[32m179\u001b[39m     result, timings = \u001b[38;5;28;43mself\u001b[39;49m\u001b[43m.\u001b[49m\u001b[43m_rpc\u001b[49m\u001b[43m(\u001b[49m\u001b[43mActionTag\u001b[49m\u001b[43m.\u001b[49m\u001b[43mEXECUTE\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mpayload\u001b[49m\u001b[43m)\u001b[49m\n\u001b[32m    180\u001b[39m \u001b[38;5;28;01mexcept\u001b[39;00m FatalError \u001b[38;5;28;01mas\u001b[39;00m e:\n\u001b[32m    181\u001b[39m     \u001b[38;5;28;01mraise\u001b[39;00m e.maybe_user_error(ir) \u001b[38;5;28;01mfrom\u001b[39;00m\u001b[38;5;250m \u001b[39m\u001b[38;5;28;01mNone\u001b[39;00m\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/hail/backend/service_backend.py:466\u001b[39m, in \u001b[36mServiceBackend._rpc\u001b[39m\u001b[34m(self, action, payload)\u001b[39m\n\u001b[32m    465\u001b[39m \u001b[38;5;28;01mdef\u001b[39;00m\u001b[38;5;250m \u001b[39m\u001b[34m_rpc\u001b[39m(\u001b[38;5;28mself\u001b[39m, action: ActionTag, payload: ActionPayload) -> Tuple[\u001b[38;5;28mbytes\u001b[39m, \u001b[38;5;28mstr\u001b[39m]:\n\u001b[32m--> \u001b[39m\u001b[32m466\u001b[39m     \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[38;5;28;43mself\u001b[39;49m\u001b[43m.\u001b[49m\u001b[43m_cancel_on_ctrl_c\u001b[49m\u001b[43m(\u001b[49m\u001b[38;5;28;43mself\u001b[39;49m\u001b[43m.\u001b[49m\u001b[43m_async_rpc\u001b[49m\u001b[43m(\u001b[49m\u001b[43maction\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mpayload\u001b[49m\u001b[43m)\u001b[49m\u001b[43m)\u001b[49m\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/hail/backend/service_backend.py:456\u001b[39m, in \u001b[36mServiceBackend._cancel_on_ctrl_c\u001b[39m\u001b[34m(self, coro)\u001b[39m\n\u001b[32m    454\u001b[39m \u001b[38;5;28;01mdef\u001b[39;00m\u001b[38;5;250m \u001b[39m\u001b[34m_cancel_on_ctrl_c\u001b[39m(\u001b[38;5;28mself\u001b[39m, coro: Awaitable[T]) -> T:\n\u001b[32m    455\u001b[39m     \u001b[38;5;28;01mtry\u001b[39;00m:\n\u001b[32m--> \u001b[39m\u001b[32m456\u001b[39m         \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[43masync_to_blocking\u001b[49m\u001b[43m(\u001b[49m\u001b[43mcoro\u001b[49m\u001b[43m)\u001b[49m\n\u001b[32m    457\u001b[39m     \u001b[38;5;28;01mexcept\u001b[39;00m \u001b[38;5;167;01mKeyboardInterrupt\u001b[39;00m:\n\u001b[32m    458\u001b[39m         \u001b[38;5;28;01mif\u001b[39;00m \u001b[38;5;28mself\u001b[39m._batch_was_submitted:\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/hailtop/utils/utils.py:181\u001b[39m, in \u001b[36masync_to_blocking\u001b[39m\u001b[34m(coro)\u001b[39m\n\u001b[32m    179\u001b[39m task = asyncio.ensure_future(coro)\n\u001b[32m    180\u001b[39m \u001b[38;5;28;01mtry\u001b[39;00m:\n\u001b[32m--> \u001b[39m\u001b[32m181\u001b[39m     \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[43mloop\u001b[49m\u001b[43m.\u001b[49m\u001b[43mrun_until_complete\u001b[49m\u001b[43m(\u001b[49m\u001b[43mtask\u001b[49m\u001b[43m)\u001b[49m\n\u001b[32m    182\u001b[39m \u001b[38;5;28;01mfinally\u001b[39;00m:\n\u001b[32m    183\u001b[39m     \u001b[38;5;28;01mif\u001b[39;00m task.done() \u001b[38;5;129;01mand\u001b[39;00m \u001b[38;5;129;01mnot\u001b[39;00m task.cancelled():\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/nest_asyncio.py:92\u001b[39m, in \u001b[36m_patch_loop.<locals>.run_until_complete\u001b[39m\u001b[34m(self, future)\u001b[39m\n\u001b[32m     90\u001b[39m     f._log_destroy_pending = \u001b[38;5;28;01mFalse\u001b[39;00m\n\u001b[32m     91\u001b[39m \u001b[38;5;28;01mwhile\u001b[39;00m \u001b[38;5;129;01mnot\u001b[39;00m f.done():\n\u001b[32m---> \u001b[39m\u001b[32m92\u001b[39m     \u001b[38;5;28;43mself\u001b[39;49m\u001b[43m.\u001b[49m\u001b[43m_run_once\u001b[49m\u001b[43m(\u001b[49m\u001b[43m)\u001b[49m\n\u001b[32m     93\u001b[39m     \u001b[38;5;28;01mif\u001b[39;00m \u001b[38;5;28mself\u001b[39m._stopping:\n\u001b[32m     94\u001b[39m         \u001b[38;5;28;01mbreak\u001b[39;00m\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/site-packages/nest_asyncio.py:115\u001b[39m, in \u001b[36m_patch_loop.<locals>._run_once\u001b[39m\u001b[34m(self)\u001b[39m\n\u001b[32m    108\u001b[39m     heappop(scheduled)\n\u001b[32m    110\u001b[39m timeout = (\n\u001b[32m    111\u001b[39m     \u001b[32m0\u001b[39m \u001b[38;5;28;01mif\u001b[39;00m ready \u001b[38;5;129;01mor\u001b[39;00m \u001b[38;5;28mself\u001b[39m._stopping\n\u001b[32m    112\u001b[39m     \u001b[38;5;28;01melse\u001b[39;00m \u001b[38;5;28mmin\u001b[39m(\u001b[38;5;28mmax\u001b[39m(\n\u001b[32m    113\u001b[39m         scheduled[\u001b[32m0\u001b[39m]._when - \u001b[38;5;28mself\u001b[39m.time(), \u001b[32m0\u001b[39m), \u001b[32m86400\u001b[39m) \u001b[38;5;28;01mif\u001b[39;00m scheduled\n\u001b[32m    114\u001b[39m     \u001b[38;5;28;01melse\u001b[39;00m \u001b[38;5;28;01mNone\u001b[39;00m)\n\u001b[32m--> \u001b[39m\u001b[32m115\u001b[39m event_list = \u001b[38;5;28;43mself\u001b[39;49m\u001b[43m.\u001b[49m\u001b[43m_selector\u001b[49m\u001b[43m.\u001b[49m\u001b[43mselect\u001b[49m\u001b[43m(\u001b[49m\u001b[43mtimeout\u001b[49m\u001b[43m)\u001b[49m\n\u001b[32m    116\u001b[39m \u001b[38;5;28mself\u001b[39m._process_events(event_list)\n\u001b[32m    118\u001b[39m end_time = \u001b[38;5;28mself\u001b[39m.time() + \u001b[38;5;28mself\u001b[39m._clock_resolution\n",
            "\u001b[36mFile \u001b[39m\u001b[32m~/miniconda3/envs/loftee/lib/python3.11/selectors.py:561\u001b[39m, in \u001b[36mKqueueSelector.select\u001b[39m\u001b[34m(self, timeout)\u001b[39m\n\u001b[32m    559\u001b[39m ready = []\n\u001b[32m    560\u001b[39m \u001b[38;5;28;01mtry\u001b[39;00m:\n\u001b[32m--> \u001b[39m\u001b[32m561\u001b[39m     kev_list = \u001b[38;5;28mself\u001b[39m._selector.control(\u001b[38;5;28;01mNone\u001b[39;00m, max_ev, timeout)\n\u001b[32m    562\u001b[39m \u001b[38;5;28;01mexcept\u001b[39;00m \u001b[38;5;167;01mInterruptedError\u001b[39;00m:\n\u001b[32m    563\u001b[39m     \u001b[38;5;28;01mreturn\u001b[39;00m ready\n",
            "\u001b[31mKeyboardInterrupt\u001b[39m: "
          ]
        }
      ],
      "source": [
        "ht = hl.read_table(\"gs://gnomad-tmp-4day/julia/gnomad.v4.1.constraint.preprocessed.for_ac_group_computation.ht\")\n",
        "agg_ht = ht.group_by(\"mutation_type\", \"downsampling\", \"ac_group\").aggregate(\n",
        "    n=hl.struct(\n",
        "        total=hl.agg.count(),\n",
        "        most_severe_consequence=hl.agg.group_by(\n",
        "            ht.most_severe_consequence.select(\"most_severe_consequence\", \"lof\"),\n",
        "            hl.agg.count()\n",
        "        ),\n",
        "        canonical_most_severe_consequence=hl.agg.group_by(\n",
        "            ht.canonical_most_severe_consequence.select(\"most_severe_consequence\", \"lof\"),\n",
        "            hl.agg.count()\n",
        "        ),\n",
        "        mane_select_most_severe_consequence=hl.agg.group_by(\n",
        "            ht.mane_select_most_severe_consequence.select(\"most_severe_consequence\", \"lof\"),\n",
        "            hl.agg.count()\n",
        "        )\n",
        "    )\n",
        ").checkpoint(\n",
        "    \"gs://gnomad-tmp-4day/julia/gnomad.v4.1.constraint.preprocessed.for_ac_group_computation.agg_ht.ht\",\n",
        "    _read_if_exists=True,\n",
        "    #overwrite=True\n",
        ")\n",
        "agg_ht.show()\n",
        "\n"
      ]
    }
  ],
  "metadata": {
    "kernelspec": {
      "display_name": "loftee",
      "language": "python",
      "name": "python3"
    },
    "language_info": {
      "codemirror_mode": {
        "name": "ipython",
        "version": 3
      },
      "file_extension": ".py",
      "mimetype": "text/x-python",
      "name": "python",
      "nbconvert_exporter": "python",
      "pygments_lexer": "ipython3",
      "version": "3.11.2"
    }
  },
  "nbformat": 4,
  "nbformat_minor": 2
}
