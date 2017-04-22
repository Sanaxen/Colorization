set ldm=..\x64\Colorization.exe

:call exec sample\名称未設定4.bmp sample\名称未設定4_marked.bmp -o 名称未設定4_Out.bmp
:goto end

:goto 2
call exec sample\example.bmp sample\example_marked.bmp -o example_Out.bmp
call exec sample\Arch.bmp sample\Arch_marked.bmp -o Arch_Out.bmp
call exec sample\Jennifer.bmp sample\Jennifer_merked.bmp -o Jennifer_Out.bmp
call exec sample\pepper.bmp sample\pepper_marked.bmp -o pepper_Out.bmp
call exec sample\リンカーン.bmp sample\リンカーン_marked.bmp -o リンカーン_Out.bmp
call exec sample\tokyo1.bmp sample\tokyo1_marked.bmp -o tokyo1_Out.bmp

call exec sample\Arch.bmp -m sample\Arch_mask.bmp -o Arch2_Out.bmp
call exec sample\Jennifer.jpg -m sample\Jennifer_mask.bmp -o Jennifer2_Out.bmp
call exec sample\pepper.bmp -m sample\pepper_mask.bmp -o pepper2_Out.bmp

:2
call exec sample\cats.bmp sample\cats_low_m.bmp -o cats_Out.bmp
call exec sample\hair.bmp sample\hair_m.bmp -o hair_Out.bmp
call exec sample\monaco.bmp sample\monaco_low_m.bmp -o monaco_Out.bmp
call exec sample\waterfall.bmp sample\waterfall_low_m.bmp -o waterfall_Out.bmp
call exec sample\yellow.bmp sample\yellow_m.bmp -o yellow_Out.bmp
:goto end

:2
call exec sample\名称未設定1.bmp sample\名称未設定1_marked.bmp -o 名称未設定1_Out.bmp
call exec sample\名称未設定2.bmp sample\名称未設定2_marked.bmp -o 名称未設定2_Out.bmp
call exec sample\名称未設定3.bmp sample\名称未設定3_marked.bmp -o 名称未設定3_Out.bmp
call exec sample\名称未設定4.bmp sample\名称未設定4_marked.bmp -o 名称未設定4_Out.bmp

:end
copy *_Out.bmp output\*.* . /v /y
del /Q *_Out.bmp

:end