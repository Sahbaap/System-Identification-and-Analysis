# System Identification and Analysis

This is the final project I did for the Linear Dynamic System class. It is related to 
a linear MIMO system with two inputs and two outputs. There are five main tasks:

## Task #1: Model identication
Find A, B, C, D matrices based on the observed impulse response so that:

## Discriminant Functions/Decision Boundaries
* **Discriminant Functions** are functions which describe a mathematical description of the optimal decision to make. That is, they divide a decision space into optimal decisions based on some empirical criteria (min loss, min risk, etc).
* Discriminant functions partition the feature space through a set of functions $\{g_1(x),...g_k(x)\}$
* The optimal decision is $\alpha(x) = \text{argmax}\{g_1(x),...g_k(x)\}$
* For Guassian distributions (see HW 1 problem 1), we have following proof:

    \textit{Proof} Let $p(x|w_i) \approx \mathcal{N}(\mu_i, \sigma^2I)$. Then
    $$
    \begin{aligned}
        g_i(x) &= \log p(x|y=i) + \log p(y=i) \\
               &= \frac{-1}{2\sigma^2}(-2\mu_i^Tx+\mu_i^T\mu_i) + \log p(y=i) \\ 
               &= \frac{1}{\sigma^2}\mu_i^Tx+\frac{-1}{2\sigma^2}\mu_i^T\mu_i+\log p(y=i) 
    \end{aligned} 
    $$

## Bayes Risk, Bayes Error
* **Bayes Risk** is the risk of the optimal decision ($\alpha$):
$$ 
    R(\alpha) = \int R(\alpha(x)|x)p(x)dx
              = \int \sum_{y=1}^k \lambda(\alpha(x)|y)p(y|x)p(x)dx
              = \int \sum_{y=1}^k \lambda(\alpha(x)|y)p(y,x)dx
$$
* We can estimate the Bayes risk (known as **empirical risk**) by summing all the risks (wrong decisions over the testing set)
$$ \text{Given } D = \{(x_i,y_i);i=1,2,...,m\} \text{   Bayes risk = } \frac{1}{m} R(\alpha) = \sum_{j=1}^{m}\lambda(\alpha(x_j)|y=j) $$
* The **Bayesian Decision Policy** minimizes loss $\pi = \text{argmin } R_\pi$ (always picks the action with the lowest risk).
     * For the case of 0/1 loss, the Bayesian decision policy will always choose the class with the highest posterior probability $p(y=i|x)$. Accordingly, the risk of this decision is $1-p(y=i|x)$