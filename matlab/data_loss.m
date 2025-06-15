function loss_percentage = data_loss(x, y)
missing         = isnan(x) | isnan(y);
loss_percentage = sum(missing)/length(missing)*100;
end
