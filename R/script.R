library( "RProtoBuf")
setwd("/users/daniel/desktop/stAR-sampler")


##### gamma plotting #####
x = seq(0,100,0.01);
dev.off()

plot(x, dgamma(x,2,rate = 0.5), type = 'l', c(0,10))

##### new proto files ##### 

samples_proto <- RProtoBuf::readProtoFiles(files = "proto/paramdata.proto")

file_path_samples <- "cmake-build-debug/samples_serialized.bin"  
binary_data_samples <- readBin(file_path_samples, "raw", file.info(file_path_samples)$size)

msg_samples <- read(sampler_data.samples, binary_data_samples)
list_rho <- as.list(msg_samples$rho)
list_phi <- as.list(msg_samples$phi)
list_sigma_w <- as.list(msg_samples$sigma_w)
list_sigma_eps <- as.list(msg_samples$sigma_eps)
list_sigma_0 <- as.list(msg_samples$sigma_0)
list_beta <- as.list(msg_samples$beta)
list_mu0 <- as.list(msg_samples$mu_0)
#list_o <- as.list(msg_samples$o)

beta <- sapply(list_beta, function(x){x$vec_value})
mu0 <- sapply(list_mu0, function(x){x$vec_value})
rho <- sapply(list_rho, function(x){x});
phi <- sapply(list_phi, function(x){x});
sig_w <- sapply(list_sigma_w, function(x){x});
sig_eps <- sapply(list_sigma_eps, function(x){x});
sig_0 <- sapply(list_sigma_0, function(x){x});
#o <- sapply(list_o, function(x){x$vec});
#o <- sapply(o, function(x){x$vec_value});



beta_true <- c(
  0.52218,
  -1.09668,
  -0.415183,
  -1.27539,
  -0.032302
)

## beta only sample  
dev.off()
par(mfrow = c(5, 1))
par(mar = c(2, 4, 2, 1))
for (i in 1:5) {
  # Create the plot
  plot(beta[i,], type = 'l', main = "", col = 'grey', ylab = bquote(beta[.(i)]))
  mean(beta[i,])
  # Add the true beta line
  lines(1:5000, rep(beta_true[i], 5000), col = "black", type = 'l', lwd = 2)  
}
n_iter= 5000;
alpha <- 0.05;
beta_low_bound <- c();
beta_high_bound <- c();
beta_means <- c();
for(i in 1:5){
beta_low_bound <- c(beta_low_bound, sort(beta[i,])[n_iter*alpha/2]) 
beta_means <- c(beta_means, mean(beta[i,]))
beta_high_bound <- c(beta_high_bound, sort(beta[i,])[n_iter*(1-alpha/2)])
}
beta_conf_int <- data.frame(beta_low_bound, beta_means, beta_high_bound)
beta_conf_int

dev.off()


###### rho



## read rho file
dev.off()
plot(rho, type = 'l', col = 'grey', ylab = "", main = expression(rho), xlab = "Iteration")
lines(1:5000, rep(0.5, 5000), type = 'l', col = 'black', lwd = 2)
plot(rho, type = 'l')
rho_low_bound <- sort(rho)[n_iter*alpha/2]
rho_low_bound
rho_high_bound <- sort(rho)[n_iter*(1- alpha/2)]
rho_high_bound 
mean(rho)
#mu_0


## read file

mu_0_true <- c(
  1.44256,
  0.540784,
  -0.00554277,
  -1.18828,
  -0.340807,
  0.815028,
  -1.37403,
  -0.486838,
  -0.272475,
  1.639
)
dev.off()
par(mfrow = c(5, 2))
par(mar = c(2, 4, 2, 1))
for (i in 1:10) {
  plot(mu0[i,], ylab = bquote(paste(mu[0], "[", .(i), "]")) ,  type = 'l', col = 'grey');
  lines(1:5000, rep(mu_0_true[i], 5000), type = 'l',col = 'black', lwd = 2);
}

mu_low_bound <- c();
mu_high_bound <- c();
mu_means <- c();

for(i in 1:5){
  mu_low_bound <- c(mu_low_bound, sort(mu0[i,])[n_iter*alpha/2]) 
  mu_means <- c(mu_means, mean(mu0[i,]))
  mu_high_bound <- c(mu_high_bound, sort(mu0[i,])[n_iter*(1-alpha/2)])
}
mu_conf_int <- data.frame(mu_low_bound, mu_means, mu_high_bound)
mu_conf_int


##### plot of variance components

## sig_eps
dev.off()
plot(sig_eps, type = 'l', col = 'grey', main =  bquote(sigma[epsilon]^2), ylab = "", xlab = "Iteration")
lines(1:5000, rep(0.1, 5000), type = 'l', col = "black", lwd = 2)
mean(sig_eps)


## sig_w
dev.off()
plot(sig_w, type = 'l',main =  bquote(sigma[w]^2), col ='grey', ylim = c(0, 1) ,xlab = "Iteration" , ylab = "")
lines(1:5000, rep(0.9, 5000), type = 'l', col = "black", lwd =2 )
plot(sig_w, type = 'l')
mean(sig_w)

## sig_0
dev.off()
plot(sig_0[1000:5000],col = 'grey' ,type = 'l', main = bquote(sigma[0]^2), xlab = "Iteration", ylab = "")
lines(1:4000, rep(1, 4000), type = 'l', col = "green")
mean(sig_0)
plot(sig_0, type = 'l')

### phi
dev.off()
plot(phi, type = 'l', col = 'grey')
mean(phi)
lower_bound_phi <- sort(phi)[n_iter*alpha/2]
lower_bound_phi
higher_bound_phi <- sort(phi)[n_iter*(1-alpha/2)]
higher_bound_phi
#### o's
o0 <- matrix(0, ncol = 5000, nrow = 10);
col = 1;
for(i in seq(11,999811 , 200)){
  o0[,col] = o[,i];
  col = col + 1;
}
dev.off()
plot(o0[1,], type = "l")  
mean(o0[1,])
plot(o0[2,], type = "l")
mean(o0[2,])
plot(o0[3,], type = "l")
mean(o0[3,])
plot(o0[4,], type = "l")
mean(o0[4,])
plot(o0[5,], type = "l")
mean(o0[5,])
plot(o0[6,], type = "l")
mean(o0[6,])
plot(o0[7,], type = "l")
mean(o0[7,])
plot(o0[8,], type = "l")
mean(o0[8,])
plot(o0[9,], type = "l")
mean(o0[9,])
plot(o0[10,], type = "l")
mean(o0[10,])

mean <- c(
  0.790282,
  2.44915,
  3.55486,
  5.77145,
  5.61658,
  5.07802,
  7.88179,
  8.32403,
  9.91118,
  10.2209
)
mean2 <- c(
  0.799816,
  2.05917,
  0.767974,
  0.722412,
  0.146886,
  -1.35127,
  2.82488,
  -0.877237,
  -2.75092,
  1.8037
)
par(mfrow = c(5, 2))
par(mar = c(2, 4, 2, 1))
for (i in 1:10) {
  plot(o0[i,1000:5000], main = paste("o10 ", i) ,  type = 'l', col = 'grey');
  lines(1:4000, rep(mean2[i], 4000), type = 'l', col = "green");
}

