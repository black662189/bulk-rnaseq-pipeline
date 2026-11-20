//表格样式
$(document).ready( function () {
    $('table#table').DataTable({
    	 "paginate": false,
    	 "paging": false,
    	 "bSort": false,
    	 //"sScrollX":'100%',
    	// "sScrollY":'500px',
    	"bAutoWidth":false,
    });
});
// 默认隐藏带hide属性的部分
$('.hide').hide();

// 按键切换图片或表格
var a=$('[name=Button]');
a.click(function(){
    getClass=$(this).attr('class');
    getId=$(this).attr('id');   
    var b=$("div."+getClass);
    // 选择对象更改class
    b.children().each(function(){
        // 判断是否为button对应对象
        if($(this).attr("id") == getId){
            // 是 则更改class为active
            b.children("#"+getId).removeClass("hide");
            b.children("#"+getId).addClass("active");
        }else{
            // 否 则更改class为hide
            b.children("#"+$(this).attr('id')).removeClass("active");
            b.children("#"+$(this).attr('id')).addClass("hide");            
        }
    $(".active").show();
    $(".hide").hide();
    });
    // console.log(b.children("#"+getId).attr('class'));
});

//sidebar收缩
var side=$('aside');
var side_width=side.width();
var section=$('section');
$('button#sidebar').click(function(){
    if(side.offset().left==0){
        side.css('margin-left',-28+'%');
        // 延迟响应
        // setTimeout(function(){section.animate({margin:"0% 2% 0% 2%"},1000);},200);
    }else{
        side.css('margin-left',0 + 'px');
        // section.animate({margin:"0% 2% 0% 25%"},1000);
    }
})

// console.log($('body').width()-side_width);
// console.log($('section').animate());
