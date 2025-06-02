set.seed(101)
N <- 10000
alpha <- runif(N, 0, 5)
beta <- runif(N, -5, 0)
x <- seq(0, 1, 0.01)
preds <- apply(cbind(alpha, beta), 1, \(beta) beta[[1]] + beta[[2]] * x) |> t()
plot(NULL, NULL, xlim = c(0, 1), ylim = c(-5, 5), xlab = "x", ylab = "y")
for (i in 1:N) {
	curve(alpha[i] + beta[i] * x, add = TRUE, col = adjustcolor("black", alpha.f = 0.05))
}

auc <- apply(preds, 1, \(y) pracma::trapz(x, y))
test <- -beta / 2
test2 <- 0.5 * 1 * alpha
