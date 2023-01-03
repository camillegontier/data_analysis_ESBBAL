using Distributions

# Matrix-based computations as explained in https://www.frontiersin.org/articles/10.3389/fncom.2016.00116/full

# Loglikelihood function for EPSP, given a set of parameters and delta_ts
function ll_binomial_tau_D(EPSP,N,p,q,sigma,tau_D,delta_t)

  X=myprobamps(EPSP,N,p,q,sigma,tau_D,delta_t)
  LL=log(X[end])
  return LL

end

# Computes the left-most vector
function P_km_given_xm(km,n,p)
  q=zeros(n+1,1)
  for xm=km:n # uses fact that xm>=km if km releases occur
    # q[xm+1]=p^km*(1-p)^(xm-km)*binomial(xm,km)
    q[xm+1]=pdf(Binomial(xm,p),km)
  end
  return q
end

# Computes the middle matrix
function P_xmp1_and_km_given_xm(km,n,p,g)
  q2=P_xmp1_given_km_and_xm(km,n,g)
  q1=P_km_given_xm(km,n,p)
  n=length(q1)-1
  for x=0:n
    q2[:,x+1]=q2[:,x+1]*q1[x+1]
  end
  return q2
end

# g is the refilling probability
function P_xmp1_given_km_and_xm(km,n,g)
  q=zeros(n+1,n+1)
  for xm=km:n
    for xmp1=xm-km:n
      # q[xmp1+1,xm+1]=g^(xmp1-xm+km)*(1-g)^(n-xmp1)*binomial(n-xm+km,xmp1-xm+km);
      q[xmp1+1,xm+1]=pdf(Binomial(n-xm+km,g),xmp1-xm+km)
    end
  end
  return q
end

function myprobamps(EPSP,N,p,q,sigma,tau_D,delta_t)

  # unpack parameters
  tD = tau_D

  A = EPSP
  ts = delta_t

  T=length(ts);
  Pprod=zeros(T,1)

  # the right vector, which is P(x1)
  r=zeros(N+1,1);
  r[end]=1;
  R=r;

  # and the left for m=1
  L1=zeros(1,N+1);

  for k1=0:N
    normalDistribution = Normal(q*k1, sigma)
    G1=pdf(normalDistribution,A[1])
    L1+=G1*P_km_given_xm(k1,N,p)';
  end
  Pprod[1]=(L1*R)[1];

  # needed for the ongoing calculation
  QR=R;
  for m=1:T-1

    g=1-exp(-delta_t[m+1]/tD);

    Lm=zeros(1,N+1);
    Qm=zeros(N+1,N+1);

    for km=0:N

      # get the Gaussian factor and matrix
      normalDistribution = Normal(q*km, sigma)
      Gm=pdf(normalDistribution,A[m]);
      Qm=Qm+Gm*P_xmp1_and_km_given_xm(km,N,p,g);

      # now the Gaussian and left vector
      Gmp1=pdf(normalDistribution,A[m+1]);
      Lm=Lm+Gmp1*P_km_given_xm(km,N,p)';

    end

    QR=Qm*QR
    Pprod[m+1]=(Lm*QR)[1]

  end

  return Pprod

end
