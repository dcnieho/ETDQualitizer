function loss_percentage = data_loss_nominal(x, y, duration, frequency)
N_valid         = sum(~(isnan(x) | isnan(y)));
loss_percentage = (1-N_valid/(duration*frequency))*100;
