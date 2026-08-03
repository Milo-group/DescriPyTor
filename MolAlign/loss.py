import torch

def loss_function (output, target, sigma = 20):

    MAE = torch.pow(torch.sum(torch.pow (target - output, 2), 1),0.5)
    MSE = torch.sum(torch.pow (target - output, 2), 1)
    SIG = 1 - 2 / (1 + torch.exp (sigma*MAE)) #20,30

    return SIG + 5*MAE #+ MSE


def loss_calc(output, target, sigma = 20):
    
    loss = 0
    
    for point in output:
        loss += torch.min(loss_function(point, target, sigma = 20))
        
    return (loss)
