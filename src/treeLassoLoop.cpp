#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <cstdio>
#include <ctime>


// [[Rcpp::export]]
double cal2norm(const arma::mat A, const arma::umat gIdx, const bool Atranspose , const int TgCB_T ) {

    double gnormsum = 0.;
    arma::mat AA;
    if(TgCB_T){
        int p = A.n_rows;
        AA = arma::zeros<arma::mat>( p*TgCB_T, A.n_cols );
        for( unsigned int i=0; i<TgCB_T; ++i )
        {
            AA.rows( arma::span(p*i, p*(i+1)-1) ) = A;
        }
    }

    for( unsigned int i=0; i<gIdx.n_rows; ++i ){
        if( !Atranspose ){
            if(!TgCB_T){
                gnormsum += std::sqrt( arma::as_scalar( arma::accu ( arma::square(A.rows(gIdx(i,0)-1, gIdx(i,1)-1)) ) )) ;
            }else{
                gnormsum += std::sqrt( arma::as_scalar( arma::accu ( arma::square(AA.rows(gIdx(i,0)-1, gIdx(i,1)-1)) ) )) ;
            }
        }else{
            gnormsum += std::sqrt( arma::as_scalar( arma::accu ( arma::square(A.cols(gIdx(i,0)-1, gIdx(i,1)-1)) ) )) ;
        }
    }

  return gnormsum;
}

// [[Rcpp::export]]
arma::mat shrink( const arma::mat A, const arma::umat gIdx, const bool Atranspose, const int TgCB_T ) {

    arma::mat R; //gnormMat;
    arma::uvec idx;
    arma::vec gnorm;
    int p = A.n_rows;
    arma::mat AA;

    if(!TgCB_T){
        R = arma::zeros<arma::mat>( p, A.n_cols );
    }else{
        R = arma::zeros<arma::mat>( p*TgCB_T, A.n_cols );
        AA = R;
        for( unsigned int i=0; i<TgCB_T; ++i )
        {
            AA.rows( arma::span(p*i, p*(i+1)-1) ) = A;
        }
    }
    
    for( unsigned int i=0; i<gIdx.n_rows; ++i )
    {
        if( !Atranspose  ){
            if(!TgCB_T){
                gnorm = arma::conv_to<arma::vec>::from( arma::sqrt( arma::sum ( arma::square(A.rows(gIdx(i,0)-1, gIdx(i,1)-1)) , 0 ) ) );
            }else{
                gnorm = arma::conv_to<arma::vec>::from( arma::sqrt( arma::sum ( arma::square(AA.rows(gIdx(i,0)-1, gIdx(i,1)-1)) , 0 ) ) );
            }
        }else{
            gnorm = arma::sqrt( arma::sum ( arma::square(A.cols(gIdx(i,0)-1, gIdx(i,1)-1)) , 1 ) );
        }

//        gnorm.elem( arma::find(gnorm < 1.) ).ones();
//        gnormMat = gnorm;
//        gnormMat.resize( gnorm.n_elem, gIdx(i,2) );
        for(unsigned int j=0; j<gIdx(i,2); ++j){
            //            gnormMat.col(j) = gnorm;
            
            if( gnorm(j) < 1. ) gnorm(j) = 1.;
            
            if( !Atranspose ){
                if(!TgCB_T){
                    R.submat( gIdx(i,0)-1, j, gIdx(i,1)-1, j ) = A.submat( gIdx(i,0)-1, j, gIdx(i,1)-1, j ) / gnorm(j);//Mat.t();
                    
//                    R( arma::span(gIdx(i,0)-1, gIdx(i,1)-1) , arma::ones<arma::uvec>(1)*j ) = A( arma::span(gIdx(i,0)-1, gIdx(i,1)-1) , arma::ones<arma::uvec>(1)*j ) / gnorm(j);//Mat.t();
                }else{
                    R.submat( gIdx(i,0)-1, j, gIdx(i,1)-1, j ) = AA.submat( gIdx(i,0)-1, j, gIdx(i,1)-1, j ) / gnorm(j);//Mat.t();
                }
            }else{
                R.submat( j, gIdx(i,0)-1, j, gIdx(i,1)-1 ) = A.submat( j, gIdx(i,0)-1, j, gIdx(i,1)-1 ) / gnorm(j); //gnormMat;
            }
        }

    }
    
  return R;
}

arma::mat TglassoCB( const arma::mat B, const int T ) {

    int p = B.n_rows/T;
    int m = B.n_cols;
    arma::mat tmp;
//    arma::vec tmp_ex;
    arma::mat CB = arma::zeros<arma::mat>( p, m );//( p*T, m );
//    arma::uvec TgIdx = arma::linspace<arma::uvec>( 0, (T-1)*p, T );
    
    for( unsigned int i=0; i<m; ++i )
    {
        /*for( unsigned int j=0; j<p; ++j )
        {
            CB(j,i) = arma::accu( B( TgIdx+j, arma::ones<arma::uvec>(1)*i ) );
        }*/
        
        tmp = B.col( i );
        tmp.reshape(p, T);
//        CB(arma::span(0,p-1), arma::span(i,i))
//        tmp_ex = arma::sum( tmp, 1);
//        tmp_ex.resize( p*T );
        CB.col(i) = arma::sum( tmp, 1);//tmp_ex;

    }
    
    /*for( unsigned int i=1; i<T; ++i )
    {
        CB.rows( arma::span(p*i, p*(i+1)-1) ) = CB.rows( arma::span(0,p-1) );
    }*/
  
  return CB;
}

// create a list of nonzero coefficients
arma::umat createGammaMask( arma::mat Beta, const int nFixedPredictors, const int nOutcomes )
{
    
    // CREATE HERE THE GAMMA "MASK"
    // INITIALISE THE INDEXES FOR THE GAMMA MASK
    arma::umat mask = arma::zeros<arma::umat>(nFixedPredictors*nOutcomes,2); //this is just an initialisation
    for( unsigned int j=0; j<nFixedPredictors; ++j)
    {
        for(unsigned int k=0 ; k<nOutcomes ; ++k)  //add gammas for the fixed variables
        {
            mask(j*nOutcomes+k,0) = j; mask(j*nOutcomes+k,1) = k;
        }
    }
    if( nFixedPredictors > 0 )
        Beta.shed_rows( 0, nFixedPredictors-1);
    for( unsigned int k=0 ; k<nOutcomes ; ++k )  //add the other gammas
    {
        arma::uvec tmpUVec = arma::find(Beta.col(k) != 0.);
        unsigned int tmpIdx = mask.n_rows;
        
        if( tmpUVec.n_elem > 0 )
        {
            mask.insert_rows( tmpIdx , arma::zeros<arma::umat>( tmpUVec.n_elem , 2 ));
            mask.submat( tmpIdx, 0, mask.n_rows-1 , 0 ) = tmpUVec + nFixedPredictors ;
            mask.submat( tmpIdx, 1, mask.n_rows-1 , 1 ).fill(k);
        }
    }
    // Gamma mask done
    
    return mask;
}

arma::mat createXB( const arma::mat& X, const arma::mat  Beta, const arma::umat&  GammaMask, const int nObservations, const int nOutcomes )
{
    arma::uvec singleIdx_k(1), VS_IN_k;
    arma::mat XB = arma::zeros<arma::mat>(nObservations,nOutcomes);
    
    if(GammaMask.n_rows > 0)
    {
        for(unsigned int k=0; k<nOutcomes; ++k)
        {
            singleIdx_k(0) = k;
            VS_IN_k =  GammaMask( arma::find(  GammaMask.col(1) == k ) , arma::zeros<arma::uvec>(1) );
            XB.col(k) = X.cols(VS_IN_k) *  Beta.submat(VS_IN_k,singleIdx_k);
        }
    }
    
    return XB;
}

arma::mat createBX( const arma::mat  Beta, const arma::mat& X, const arma::umat&  GammaMask )
{
    arma::uvec singleIdx_k(1), VS_IN_k;
    arma::mat BX = arma::zeros<arma::mat>(Beta.n_rows, X.n_cols);
    
    if(GammaMask.n_rows > 0)
    {
        arma::uvec VS_IN_tmp = arma::unique(GammaMask.col(0));
//        Rcpp::Rcout << "VS_IN_tmp=" << VS_IN_tmp.t();
        for(unsigned int k=0; k<VS_IN_tmp.n_elem; ++k)
        {
            singleIdx_k(0) = k;
            VS_IN_k =  GammaMask( arma::find(  GammaMask.col(0) == VS_IN_tmp(k) ) , arma::ones<arma::uvec>(1) );
//            Rcpp::Rcout << "size(B)=" << arma::size(Beta) << "; size(C)=" << arma::size(X) << "; k=" << k;
            BX.row(k) = Beta.submat(singleIdx_k,VS_IN_k) * X.rows(VS_IN_k);
        }
    }
    
    return BX;
}

// [[Rcpp::export]]
arma::vec XXeigen( const arma::mat& X, const arma::mat& V, const arma::uvec& tIdx ) {
    
    int T = arma::max( tIdx );
    arma::vec L = arma::zeros<arma::vec>( T );
    arma::uvec tIdx_i;
    arma::mat XX_i;
    
    for( unsigned i=0; i<T; ++i ){
            
        tIdx_i = arma::find( tIdx == i+1 );
        XX_i = X.rows( tIdx_i ).t() * V( tIdx_i, tIdx_i ) * X.rows( tIdx_i );

        arma::sp_mat tmp( XX_i );
        XX_i.clear();
            
        L(i) = arma::as_scalar( arma::real(arma::eigs_gen( tmp, 1 )) );
//        Rcpp::Rcout << "L(" << i+1 << ")=" << L(i) << "\n";
            
        Rcpp::checkUserInterrupt(); // this checks for interrupts from R
    }
    
    return L;
}

    
// [[Rcpp::export]]
arma::sp_mat mixLassoLoop(const arma::mat& V, const arma::sp_mat& C, const arma::umat& gIdx, const bool Tglasso, const arma::uvec& tIdx, const arma::mat& X, const arma::mat& Y, const bool intercept, const int num_nonpen, arma::vec L, const double lambda, const arma::vec& option, const double mu, const int NoVar, const double gamma, const arma::umat& y_mis, const double alpha ) {

    double obj = 0.;
    double obj_new = 0.;
    
    
    double theta = 1.;
    double theta_new = 1.;
    int nOutcomes = Y.n_cols;
//    int nObservations = Y.n_rows;
    int nFixedPredictors = num_nonpen + int(intercept);
    int p = X.n_cols;
    int T = arma::max( tIdx );
    arma::uvec predictorFirst(T), predictorLast(T), tIdx_i;

    arma::vec num_nonzero = arma::zeros<arma::vec>(NoVar); // monitor the number of nonzero coefficients over the last 10 iterations
    arma::umat GammaMask;
    arma::mat XB;
    //arma::sp_mat XX_i, XY_i;
    
    arma::mat RC_tree, CBt, TgCBt;
//    arma::sp_mat /*Tglasso_C,*/ RC_Tglasso;
    arma::mat RC_Tglasso;
    arma::umat Tglasso_gIdx;
    
    arma::rowvec bw0, bv0;
    arma::mat bW0, bw1, bv, bx_update;
    arma::mat grad_bw0, grad_bW0, grad_bw;
    arma::mat x, res;
    arma::uvec vsIdx;
    arma::uvec bw1_nonzero;
    //arma::mat bx = arma::zeros<arma::mat>(bx.n_rows, bx.n_cols);
    arma::mat bx;
    
    Rcpp::Rcout <<"start...\n";
//    std::clock_t start;
//    double duration1, duration2, duration3, duration4, duration5;
//    start = std::clock();
    
    arma::cube XX_i(p, p, T), XY_i(p, Y.n_cols, T);
//    arma::mat XX_i, XY_i;
//    arma::vec eigval;
    for( unsigned i=0; i<T; ++i ){
        
        tIdx_i = arma::find(tIdx==i+1);
        XX_i.slice(i) = X.rows( tIdx_i ).t() * V( tIdx_i, tIdx_i ) * X.rows( tIdx_i );
        XY_i.slice(i) = X.rows( tIdx_i ).t() * V( tIdx_i, tIdx_i ) * Y.rows( tIdx_i );
        
        if( L(0)==0. ){
            arma::sp_mat tmp( XX_i.slice(i) );
            //        arma::sp_mat tmp( XX_i );
            //        Rcpp::Rcout << "XX_i: \n" << XX_i.subcube(0,0,i,4,4,i) <<  "\n";
                    
            L(i) = arma::as_scalar( arma::real(arma::eigs_gen( tmp, 1 )) ) + L(i);
        }
        
        Rcpp::checkUserInterrupt(); // this checks for interrupts from R
    }
//    L0 = L;
    
//    Rcpp::Rcout << "L(i)=" << L.t() <<  "\n";
    if( intercept ){
        vsIdx = arma::linspace<arma::uvec>(0, (p+1)*T-1, (p+1)*T);
        vsIdx.shed_rows( arma::linspace<arma::uvec>(0, (p+1)*T-1, T) );
        bx = arma::zeros<arma::mat>((1+p)*T, Y.n_cols);
    }else{
        vsIdx = arma::linspace<arma::uvec>(0, p*T-1, p*T);
        vsIdx.shed_rows( arma::linspace<arma::uvec>(0, p*T-1, T) );
        bx = arma::zeros<arma::mat>(p*T, Y.n_cols);
    }
    
    if( C.n_rows > vsIdx.n_elem ){
        CBt = C * (bx.rows(vsIdx)).t();
        RC_tree = (shrink(CBt, gIdx, false, 0)).t() * C;
    }else{
        CBt = bx.rows(vsIdx) * C.t();
        RC_tree = shrink(CBt, gIdx, true, 0) * C;
    }
    
    if( Tglasso ){
        // construct gIdx for the grouped tissue features
        Tglasso_gIdx = arma::ones<arma::umat>( p-num_nonpen, 3 ) * T;
        Tglasso_gIdx.col(1) = arma::linspace<arma::uvec>( 1, p-num_nonpen, p-num_nonpen ) * T;
        Tglasso_gIdx.col(0) = Tglasso_gIdx.col(1) - (T-1);

        //TgCBt = Tglasso_C * bx.rows(vsIdx) / mu;
        TgCBt = TglassoCB( bx.rows(vsIdx), T ) * (1.-alpha) * gamma * std::sqrt(T)/ mu; /* use a trick to avoid Tglasso_C which is a too big matrix */
        //RC_Tglasso = Tglasso_C.t() * shrink(TgCBt, Tglasso_gIdx, false);
        RC_Tglasso = TglassoCB( shrink(TgCBt, Tglasso_gIdx, false, T), T );
        
    }else{
        RC_Tglasso = 0.;
    }
    
    for( unsigned int i=0; i<T; ++i ){

        if( intercept ){
            predictorFirst(i) = i*(1+p);
            predictorLast(i) = i*(1+p)+p;
        }else{
            predictorFirst(i) = i*p;
            predictorLast(i) = i*p+p-1;
        }
        
    }
    
    for( unsigned int iter=0; iter<option(0); ++iter ){

        Rcpp::checkUserInterrupt(); // this checks for interrupts from R
        
        theta_new = 2./(double(iter)+3.);
        obj_new = 0.;
        //obj_approx = 0.;
        
        // update tissue-specific parameters
        for( unsigned int i=0; i<T; ++i ){

            tIdx_i = arma::find(tIdx==i+1);
            bw1 = bx.rows(predictorFirst(i), predictorLast(i));

            // run onegrad for each tissue information
            if( intercept ){
                bw0 = bw1.row(0);
                if( num_nonpen == 0 ){
                    bw1.shed_row(0);
//                    arma::sp_mat tmp(bw1);

                    GammaMask = createGammaMask( bw1, 0, nOutcomes );
                    
//                    grad_bw0 = arma::ones<arma::rowvec>( tIdx_i.n_elem ) * ( arma::ones<arma::vec>( tIdx_i.n_elem) * bw0 + X.rows( tIdx_i ) * tmp - Y.rows( tIdx_i ) );
                    grad_bw0 = arma::ones<arma::rowvec>( tIdx_i.n_elem ) * ( arma::ones<arma::vec>( tIdx_i.n_elem) * bw0 + createXB( X.rows( tIdx_i ), bw1, GammaMask, tIdx_i.n_elem, nOutcomes ) - Y.rows( tIdx_i ) );

                    if( !Tglasso ){
//                        grad_bw = XX_i.slice(i) * tmp - XY_i.slice(i) + RC_tree.rows(i*p, i*p+p-1);
                        grad_bw = createXB( XX_i.slice(i), bw1, GammaMask, tIdx_i.n_elem, nOutcomes ) - XY_i.slice(i) + RC_tree.rows(i*p, i*p+p-1);
                    }else{
                        
//                        grad_bw = XX_i.slice(i) * tmp - XY_i.slice(i) + RC_tree.rows(i*p, i*p+p-1) + RC_Tglasso.rows(i*(p-num_nonpen),i*(p-num_nonpen)+p-num_nonpen-1);
//                        XB = createXB( XX_i.slice(i), bw1, GammaMask, p, nOutcomes );
                        grad_bw = createXB( XX_i.slice(i), bw1, GammaMask, p, nOutcomes ) - XY_i.slice(i) + RC_tree.rows(i*p, i*p+p-1) + RC_Tglasso;//.rows(i*(p-num_nonpen),i*(p-num_nonpen)+p-num_nonpen-1);
                    }
                }else{
                    bW0 = bw1.rows(1,num_nonpen);
                    bw1.shed_rows(0, num_nonpen);
                    
                    grad_bw0 = arma::ones<arma::rowvec>( tIdx_i.n_elem ) * ( arma::ones<arma::vec>( tIdx_i.n_elem ) * bw0 + X.rows( tIdx_i ) * join_cols(bW0,bw1) - Y.rows( tIdx_i ) );
                    grad_bW0 = (X.cols(0,num_nonpen-1)).t() * ( arma::ones<arma::vec>( tIdx_i.n_elem )*bw0 + X.rows( tIdx_i )*join_cols(bW0,bw1) - Y.rows( tIdx_i ));
                    
                    x = X.rows( tIdx_i );
                    x.shed_cols(0,num_nonpen-1);
                    grad_bw = x.t() * ( arma::ones<arma::vec>( tIdx_i.n_elem )*bw0 + X.rows( tIdx_i )*join_cols(bW0,bw1) - Y.rows( tIdx_i )) + RC_tree + RC_Tglasso;
                }
                // update regression intercepts
                bw0 = bw0 - 1./L(i) * grad_bw0;
            }else{
                if( num_nonpen > 0 ){
                    bW0 = bx.rows(0,num_nonpen-1);
                    bw1.shed_rows(0, num_nonpen-1);
                }
            }

            if( num_nonpen > 0 )
                bW0 = bW0 - grad_bW0/L(i);
            
            bv = bw1 -  grad_bw/L(i);
            //bw1 = abs( bv ) - (1.+l1_alpha*alpha) * lambda/L(i); /* do l1-penalty more when adding group tissue penalty*/
            //bw1 = abs( bv ) - lambda/L(i);
            bw1 = abs( bv ) - (alpha*gamma + lambda)/L(i);
            bw1.elem( arma::find(bw1<0) ).zeros();
            bw1 = sign( bv ) % bw1;
            
//            obj_approx += arma::accu( arma::square(bv0 - bw0) ) + arma::accu( arma::square(bw1 - bv) ) + (alpha*gamma + lambda)/L(i) * arma::accu(arma::abs(bw1));
            
            if( num_nonpen == 0 ){
                bx_update = join_cols( bw0, bw1 );
            }else{
                bx_update = join_cols( bw0, bW0, bw1 );
            }
            
            GammaMask = createGammaMask( bx_update, nFixedPredictors, nOutcomes );
            
            bv.clear();
            bw1.clear();
            
            // calculate residuals
            if( intercept ){
//                res = Y.rows( tIdx_i ) - join_rows(arma::ones<arma::vec>(tIdx_i.n_elem),X.rows( tIdx_i )) * bx_update;
                res = Y.rows( tIdx_i ) - createXB( join_rows(arma::ones<arma::vec>(tIdx_i.n_elem),X.rows( tIdx_i )), bx_update, GammaMask, tIdx_i.n_elem, nOutcomes );
                
                if( arma::accu(y_mis) > 0 ){
                     obj_new += arma::accu( square( res.elem( arma::find(y_mis.rows( tIdx_i )==0) ) ) )/2./tIdx_i.n_elem;// /Y.n_cols
                }else{
                    res = res.t() * V( tIdx_i, tIdx_i ) * res;
                    obj_new += arma::accu( res.diag() )/2./tIdx_i.n_elem;
                }
//                    Rcpp::Rcout << "res" << i << "=" << arma::accu( square( res)) << "\n";
                
            }else{
                res = Y.rows( tIdx_i ) - X.rows( tIdx_i ) * bx_update;
                obj_new += arma::accu( square( res.elem( arma::find(y_mis.rows( tIdx_i )==0) ) ) )/2./tIdx_i.n_elem;// /Y.n_cols
            }
            
            bx.rows(predictorFirst(i), predictorLast(i)) = bx_update + (1.-theta)/theta * theta_new * ( bx_update-bx.rows(predictorFirst(i), predictorLast(i)) );
            
        }
        res /= T;
        
        if( C.n_rows > vsIdx.n_elem ){
            CBt = C * (bx.rows(vsIdx)).t();
            RC_tree = (shrink(CBt / mu, gIdx, false, 0)).t() * C;
            obj_new += cal2norm(CBt, gIdx, false, 0);
        }else{
            CBt = bx.rows(vsIdx)*C.t();
            RC_tree = shrink(CBt / mu, gIdx, true, 0) * C;
            
            obj_new += cal2norm(CBt, gIdx, true, 0);
        }

        if( Tglasso ){
            //TgCBt = Tglasso_C*bx.rows(vsIdx) / mu;
            TgCBt = TglassoCB( bx.rows(vsIdx), T ) * (1.-alpha) * gamma * std::sqrt(T)/ mu;
            obj_new += cal2norm(TgCBt, Tglasso_gIdx, false, T);
            
            //RC_Tglasso = Tglasso_C.t() * shrink(TgCBt, Tglasso_gIdx, false);
            RC_Tglasso = TglassoCB( shrink(TgCBt, Tglasso_gIdx, false, T), T );
        }
        
        // print something on how the iteration is going
        // if( (iter+1) % 30 == 0 )
            Rcpp::Rcout << "iter=" << iter+1 << "; (obj_new-obj)/obj=" << std::abs((obj_new-obj)/obj) << "; res^2=" << arma::accu( res.diag() ) << "\n";
        
        if( (iter>10) && (std::abs((obj_new-obj)/obj)<option(2)) )
            break;
        theta = theta_new;
        obj = obj_new;
        
        bw1_nonzero = arma::find( bx!=0 );//arma::abs(bx)>0.0001 );
        Rcpp::Rcout <<"bw1_nonzero: "<< bw1_nonzero.n_elem //<<"; sd(num.nonzero)="<< arma::stddev(num_nonzero)
        << "; max(b)=" << (arma::abs(bx)).max() << "; median(b)=" << arma::median(arma::median(arma::abs(bx(bw1_nonzero)))) << '\n';
        
        // monitor the number of nonzero coefficients over the last 10 iterations
        if( iter < NoVar ){
            num_nonzero(iter) = bw1_nonzero.n_elem;
        }else{
            num_nonzero.subvec(0,NoVar-2) = num_nonzero.subvec(1,NoVar-1);
            num_nonzero(NoVar-1) = bw1_nonzero.n_elem;
            
            if( arma::stddev(num_nonzero) == 0. ){
                Rcpp::Rcout << "Stopped due to no variation of the number of nonzero coefficients over the last " << NoVar << " iterations!\n";
                break;
            }
        }
        
        // stop if the estimated model is too dense or no variation of number of nonzero bx over the last 10 iterations
        if( /*((p*T*Y.n_cols>1000000) && (bw1_nonzero.n_elem > 0.1*p*T*Y.n_cols) && (iter > 30)) | ((iter>100) && (arma::stddev(num_nonzero)==0.)) |*/
           ((bw1_nonzero.n_elem > 0.5*p*T*Y.n_cols) && (iter > 30)) |
           (bw1_nonzero.n_elem>1000000) |
           ((arma::abs(bx)).max()>50) ){
            bx.fill(999.);
            Rcpp::Rcout << "Stopped because the estimated model is too dense or large betas!\n";
            break;
        }
        
        if( iter == option(0) )
            Rcpp::Rcout << "Stopped due to the setting of maximum iterations!\n";

    }
    
    
    return arma::sp_mat(bx);
    //return Rcpp::List::create(Rcpp::Named("bx") = arma::sp_mat(bx),
    //                        Rcpp::Named("obj_new") = obj_new);
}

// [[Rcpp::export]]
arma::sp_mat treeLassoLoop(const arma::mat& X, const arma::mat& Y, const arma::sp_mat& C, const arma::umat& gIdx, const double TauNorm, const bool intercept, const int num_nonpen, const double lambda, const arma::vec& option, const double mu) {
    
    double obj = 0.;
    double obj_new = 0.;
    double L;
    arma::mat XX, XY, R, CB, res;
    
    XX = X.t() * X;
    XY = X.t() * Y;
    
    arma::sp_mat tmp(XX);
    L = arma::as_scalar( arma::real(arma::eigs_gen( tmp, 1 )) ) + lambda*lambda * TauNorm/mu;
    
    int nOutcomes = Y.n_cols;
    int nObservations = Y.n_rows;
    int p = X.n_cols;
    
    arma::rowvec bw0;
    arma::mat bW0, bw1, bv, bx_update;
    arma::mat grad_bw0, grad_bW0, grad_bw;
    arma::mat bx;
    
    // initialize bw0, b0, bw1, bW0, bx
    if( intercept )
        bw0 = arma::zeros<arma::rowvec>(nOutcomes);

    if( num_nonpen == 0 ){
        bw1 = arma::zeros<arma::mat>(p, nOutcomes);
        if( intercept ){
            bx = join_cols( bw0, bw1 );
        }else{
            bx = bw1;
        }
    }else{
        bW0 = arma::zeros<arma::mat>(num_nonpen, nOutcomes);
        bw1 = arma::zeros<arma::mat>(p - num_nonpen, nOutcomes);
        if( intercept ){
            bx = join_cols( bw0, bW0, bw1 );
        }else{
            bx = join_cols( bW0, bw1 );
        }
    }
    
    Rcpp::Rcout <<"start...\n";

    double theta = 1.;
    double theta_new = 1.;
    
    
    for( unsigned int iter=0; iter<option(0); ++iter ){

        Rcpp::checkUserInterrupt(); // this checks for interrupts from R
        
        theta_new = 2./(double(iter)+2.);
        obj_new = 0.;
        
        // compute grad(f(w_k))
        R = shrink( C * bw1.t()/mu, gIdx, false, 0);
    
        if( num_nonpen==0 ){
            if( intercept){
                grad_bw0 = arma::zeros<arma::rowvec>(nObservations) * (arma::zeros<arma::vec>(nObservations)*bw0 + X*bw1 - Y);
            }
            grad_bw = XX * bw1 - XY + R.t() * C;
        }else{
            if( intercept ){
                
                grad_bw0 = arma::zeros<arma::rowvec>(nObservations) * (arma::zeros<arma::vec>(nObservations)*bw0 + X*join_cols(bW0,bw1) - Y);
            }
            
            grad_bW0 = X.cols(0,num_nonpen-1).t() * (arma::zeros<arma::vec>(nObservations) *bw0 + X*join_cols(bW0,bw1) - Y);
            grad_bw = X.cols(num_nonpen,p-1).t() * (arma::zeros<arma::vec>(nObservations) *bw0 + X*join_cols(bW0,bw1) - Y) + R.t() * C;
            
        }
        
        bv = bw1 -  grad_bw/L;
        bw1 = abs( bv ) - lambda/L;
        bw1.elem( arma::find(bw1<0) ).zeros();
        bw1 = sign( bv ) % bw1;
        
        if( num_nonpen == 0 ){
            if( intercept ){
                bw0 = bw0 - grad_bw0/L;
                bx_update = join_cols( bw0, bw1 );
            }else{
                bx_update = bw1;
            }
        }else{
            bW0 = bW0 - grad_bW0/L;
            if( intercept ){
                bw0 = bw0 - grad_bw0/L;
                bx_update = join_cols( bw0, bW0, bw1 );
            }else{
                bx_update = join_cols( bW0, bw1 );
            }
        }
        
        //bv.clear();
        //bw1.clear();
        
        // calculate residuals
        if( intercept ){
            res = Y - join_rows(arma::ones<arma::vec>(nObservations),X) * bx_update;
        }else{
            res = Y - X * bx_update;
        }
        
        CB = C * bx_update.t();
        obj_new = arma::accu( square( res/*.elem( arma::find(y_mis==0) )*/ ) )/2./nObservations + cal2norm(CB, gIdx, false, 0);
        
        bx = bx_update + (1.-theta)/theta * theta_new * ( bx_update-bx);
        
        if( num_nonpen>0 ){
            if( intercept ){
                bw0 = bx.row(0);
                bW0 = bx.rows(1,num_nonpen);
            }else{
                bW0 = bx.rows(0,num_nonpen-1);
                bw1 = bx.rows(num_nonpen, p-1);
            }
        }else{
            if(intercept){
                bw0 = bx.row(0);
                bw1 = bx.rows(1,p);
            }else{
                bw1 = bx;
            }
        }
        
        // print something on how the iteration is going
         if( (iter+1) % 50 == 0 )
            Rcpp::Rcout << "iter=" << iter+1 << "; (obj_new-obj)/obj=" << std::abs((obj_new-obj)/obj) << "\n";
        

        if( (iter>10) && (std::abs((obj_new-obj)/obj)<option(2)) )
            break;
        theta = theta_new;
        obj = obj_new;
        
        
        if( iter == option(0) )
            Rcpp::Rcout << "Stopped due to the setting of maximum iterations!\n";

    }
    
    
    return arma::sp_mat(bx);
    //return Rcpp::List::create(Rcpp::Named("bx") = arma::sp_mat(bx),
    //                        Rcpp::Named("obj_new") = obj_new);
}
/*
double pre.grad(const arma::umat Tree, const arma::vec Tw, const bool Tglasso ) {

    arma::umat gIdx;
    
    int V = Tree.n_rows;
    int K = Tree.n_cols;
    
    int SV = arma::accu( Tree );
    arma::uvec sum_col_T = arma::sum( Tree, 1 );
    arma::uvec csum = cumsum( sum_col_T );
    
    if( csum.n_elem != 1 ){
        gIdx = join_rows( join_cols(1, csum.subvec(0,csum.n_elem-2)+1), csum, sum_col_T );
    }else{
        gIdx = join_cols(1, csum, sum_col_T );
    }
    
    arma::uvec J = arma::zeros<arma::uvec>( SV );
    arma::uvec W = J;
    
    for( unsigned int v=0; v<V-1; v++ ){
        J.subvec( gIdx(v,0), gIdx(v,1) ) = arma::find( Tree.row(v)==1 );
        W.subvec( gIdx(v,0), gIdx(v,1) ) = Tw(v);
        
        C.submat( 0, J(gIdx(v,0)), SV-1, J(gIdx(v,1)) ) =
    }
    
    arma::mat C( SV, K );
    C = arma::sp_mat<arma>(n_rows, n_cols)
    
  return Tglasso_C;
}
*/
