# Colorization
  
Visual studio 2012 C++ build
  
白黒画像のカラー化

<img src="https://github.com/Sanaxen/Colorization/blob/master/image1.png"/>
<img src="https://github.com/Sanaxen/Colorization/blob/master/image2.png"/>

共役勾配法(Conjugate Gradient Method)で収束計算しています。
精度は0.0001を保障していますが非常に遅いです。
マルチスレッドで処理速度を改善していますがSparse行列の扱いが非常に遅いです。(Colorization_old.exe)

処理速度が約30倍高速になりました。（処理ループの方法を最適化）


# 比較用のビルドに関して
  
Anat_Levin版
  
#define _used_Anat_Levin_VERSION

注）<http://www.cs.huji.ac.il/~yweiss/Colorization/>
著作権：2004年のヘブライ大学エルサレム大学。全著作権所有
マルチグリッド法で計算していますがマルチグリッド法部分は研究用に無償公開されているコードをリンクしています。

>This package contains an implementation of the image colorization approach described in the paper:
>A. Levin D. Lischinski and Y. Weiss Colorization using Optimization.
>ACM Transactions on Graphics, Aug 2004. 
> 
>
>Usage of this code is free for research purposes only. 
>Please refer to the above publication if you use the program.
>
>Copyrights: The Hebrew University of Jerusalem, 2004.
>All rights reserved.

>Written by Anat Levin.
>Please address comments/suggestions/bugs to: <alevin@cs.huji.ac.il>

