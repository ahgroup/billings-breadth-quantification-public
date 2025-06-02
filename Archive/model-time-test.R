# Time estimation

subsample_model_metadata <- targets::tar_read(subsample_model_metadata)
smml <- targets::tar_read(subsample_model_metadata_list)

source(here::here("R", "Model-Setup.R"))

inds <- with(
	subsample_model_metadata,
	which(model == "full" & censoring == "yes")
)

test_in <- smml[inds[[1]]]

brms_control_arguments <- generate_brms_control_arguments(FALSE)

test_out <- model_dispatch(test_in, brms_control_arguments)

times <- test_out$fit |> rstan::get_elapsed_time()
