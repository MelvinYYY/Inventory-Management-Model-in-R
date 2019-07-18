### (Q,R) Model for Inventory Management




QR_alpha <- function(d, sdd, lt, sdlt, u = 0, sdu = 0, k, h, c, alpha = 0.95){
  
  # Input:
  # d:    Expected Annual Demand 
  # sdd:  Standard Deviation of Demand 
  # lt:   Lead time
  # sdlt: Standard Deviation of Lead Time
  # k:    Set up Cost
  # h:    Holding cost percentage 
  # c:    Unit cost 

  # Output: a matirx of recommended order quantity, Reorder Point, Safety Stock, F(R) Standard Normal Probability, Z-score 
  
  # Main:
  
  # Calculate Demand during lead time
  u = ifelse(u == 0, 
             d * (lt/365),
             u)
  
  # Calculate Standard Deviation Demand during lead time
  sdu = ifelse(sdu == 0,
               sqrt((lt/365)*(sdd^2) + (d^2)*((sdlt/365)^2)),
               sdu)
  
  # Initial Q0 from the EOQ model
  q = round(sqrt(2*k*d/(h*c)), 0)
  
  # Find a z that satisfies loss function (z) = alpha
  z = qnorm(alpha)
  std.loss = exp(-z^2/2)/sqrt(2*pi) - z*pnorm(-z) # Standard Loss Function
  
  # Calculate inventory level 
  ss = round(z * sdu, 0) # Safety Stock
  r = round(ss + u, 0)   # Reorder Point
  
  out <- cbind(q, r, ss, alpha, z)
  
  return(out)
  
}



QR_beta <- function(d, sdd, lt, sdlt, u = 0, sdu = 0, k, h, c, p){
  
  # Input:
  # d:    Expected Annual Demand 
  # sdd:  Standard Deviation of Demand 
  # lt:   Lead time
  # sdlt: Standard Deviation of Lead Time
  # k:    Set up Cost
  # h:    Holding cost percentage 
  # c:    Unit cost 
  # p:    Unit Profit
  
  
  # Output: a matirx of recommended order quantity, Reorder Point, Safety Stock, F(R) Standard Normal Probability, Z-score 
  
  # Main:
  
  # Calculate Demand during lead time
  u = ifelse(u == 0, 
             d * (lt/365),
             u)
  
  # Calculate Standard Deviation Demand during lead time
  sdu = ifelse(sdu == 0,
               sqrt((lt/365)*(sdd^2) + (d^2)*((sdlt/365)^2)),
               sdu)
  
  # Initial Q0 from the EOQ model
  q0 = round(sqrt(2*k*d/(h*c)), 0)
  q_vec = c(q0)  # Set up vector to store value
  ss_vec = c(0)
  r_vec = c(0)   # Set up vector to store value
  fr_vec = c(0)
  z_vec = c(0)
  
  # Start irretation
  for(i in 2:30){
    
    # Calculate F(R)
    fr = 1 - (q_vec[(i-1)] * (h*c) / (p*d))
    fr_vec[i] = fr 
    z = qnorm(fr)
    z_vec[i] = z
    std.loss = exp(-z^2/2)/sqrt(2*pi) - z*pnorm(-z) # Standard Loss Function
    
    ss = round(z * sdu, 0) # Safety Stock
    ss_vec[i] = ss         # Append Data
    r = round(ss + u, 0)   # Reorder Point
    r_vec[i] = r           # Append Data
    q = round(sqrt(2*(k + (p * sdu * std.loss)) * d/(h*c)), 0) # Recommended order Quantity
    q_vec[i] = q           # Append Data
    
    # Early Stop
    if(q_vec[i] == q_vec[i-1] && r_vec[i] == r_vec[i-1]){
      break
    }
    
  }
  
  out <- cbind(q_vec, r_vec, ss_vec, fr_vec, z_vec)
  
  return(out)
  
}



