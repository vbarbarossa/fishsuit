##%######################################################%##
#                                                          #
####            Variable importance function            ####
#                                                          #
##%######################################################%##

#' Variable importance
#'
#' Calculate contribution of predictor variables to the model.
#' Function will make a reference prediction of the model using the standard set of variables.
#' Then, the values in predictor variables are randomized, and the prediction is repeated with the set of variables
#' that contain a randomized variable. Correlation coefficient is calculated between the reference prediction and randomized prediction.
#' Given importance value is \code{1 - correlation ** 2} for each variable. Number of randomizations can be set (default is one)
#'
#' @param data Input data with variables for which to calculate the variable importance. With this data you should be able to run predict function on the model.
#' @param model Model to be used for prediction. Function is tested only on glm object class.
#' @param iterations_num Number of randomization iterations. Default is 1 iteration.
#' @param clean Return cleaned data (default is \code{FALSE}). A dataframe will be returned, only with variables that participated in the model (in case of model selection).
#'
#' @return Output is a matrix where rows have variable importance value for each variable, and the columns are individual iterations. If clean = TRUE, return class is dataframe.
#' @export
#'
#' @author Mirza Cengic
#' @examples var_importance(data = mydat, model = my_model, iterations_num = 10)
#' @importFrom magrittr "%>%"
#' @importFrom tibble rownames_to_column
#' @import dplyr

variable_importance <- function(data, model, iterations_num = 1,
                                clean = FALSE)
{
  # Pass here the model and the data. Here we want to check if
  # the predictions can be calculated on the data, since the goal
  # of the function is to use the correlation between the predictor
  # and a randomized value to calculate variable importance.

  reference_prediction <- try(predict(model, data))

  if (inherits(reference_prediction, "try-error"))
  {
    stop("Error with reference prediction")
  }

  # Create matrix in which to store the values for the variable importance
  output_matrix <- matrix(0, nrow = length(names(data)),
                          ncol = iterations_num,
                          dimnames = list(names(data), paste0("Iter_", 1:iterations_num)))


  #### Loop that works (but might not be correct)
  for (iter in 1:iterations_num)
  {
    for(var_name in names(data))
    {
      # Copy the data so each iteration is independent
      dat <- data
      # print(var_name)
      # Randomize the predictor variable
      dat[, var_name] <- sample(dat[, var_name])
      # Predict on the dataset with randomized variable
      randomized_prediction <- predict(model, dat)
      # Calculate correlation between the reference and randomized prediction, and substract it from 1
      output_matrix[var_name, iter] <- 1 - round(cor(x = reference_prediction,
                                                     y = randomized_prediction,
                                                     use = "pairwise.complete.obs",
                                                     method = "pearson"), 4)
    }
  }

  if (clean)
  {
    # Get names of variables that were used for the model
    var_names <- names(model$model)
    var_names <- var_names[!var_names %in% c("PA", "(weights)")]

    output_matrix <- output_matrix %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Variable") %>%
      dplyr::filter(Variable %in% var_names)

    return(output_matrix)

  } else {
      return(output_matrix)
  }
}
